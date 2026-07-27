"""Convert coarse-grained mBuild systems to CGsmiles graphs and strings.

A coarse-grained system is either an ``mbuild.path.Path`` or an
``mbuild.Compound`` whose leaf particles are the CG beads. Both carry
the same information — bead names, bead coordinates, and a bond graph —
and both convert to the same CGsmiles meta-graph representation used by
the rest of this package.
"""

import networkx as nx
import numpy as np

__all__ = ["to_cgsmiles_graph", "to_cgsmiles"]


def to_cgsmiles_graph(system, fragname_map=None):
    """Convert a coarse-grained system to a CGsmiles-compatible meta graph.

    Each node in the returned graph carries a ``fragname`` attribute
    taken from the bead names (after applying ``fragname_map``), a
    ``beadname`` attribute holding the original bead name, and a
    ``position`` attribute. All edges are given ``order=1``.

    Parameters
    ----------
    system : mbuild.path.Path or mbuild.Compound, required
        The coarse-grained system to convert. For a ``Compound``, each
        leaf particle is treated as one CG bead: particle names become
        fragment names, particle positions become bead positions, and
        the compound's bonds become the CG bond graph.
    fragname_map : dict[str, str], optional
        Mapping of bead names to CGsmiles fragment names. Bead names not
        present in the map are used as fragment names directly. Useful
        when bead names differ from the fragment names used in the
        CGsmiles fragment string.

    Returns
    -------
    networkx.Graph
        Meta graph usable with ``cgsmiles.MoleculeResolver``.
    """
    from mbuild.compound import Compound

    if isinstance(system, Compound):
        return _compound_meta_graph(system, fragname_map)
    if hasattr(system, "bond_graph") and hasattr(system, "beads"):
        return _path_meta_graph(system, fragname_map)
    raise TypeError(
        f"Cannot convert {type(system).__name__} to a CGsmiles graph; "
        "pass an mbuild.path.Path or an mbuild.Compound whose leaf "
        "particles are the CG beads."
    )


def _fragname(name, fragname_map):
    if fragname_map is not None:
        return fragname_map.get(name, name)
    return name


def _path_meta_graph(path, fragname_map):
    """Meta graph from a Path's beads, coordinates, and bond graph."""
    # TODO: If we can set bond order in Path, grab it and set it here
    if len(path.bond_graph) != len(path.coordinates):
        path._extend_bond_graph()
    meta_graph = nx.Graph()
    for node in sorted(path.bond_graph.nodes):
        name = str(path.beads[node])
        meta_graph.add_node(
            int(node),
            fragname=_fragname(name, fragname_map),
            beadname=name,
            position=np.asarray(path.coordinates[node], dtype=float),
        )
    for u, v in path.bond_graph.edges:
        meta_graph.add_edge(int(u), int(v), order=1)
    return meta_graph


def _compound_meta_graph(compound, fragname_map):
    """Meta graph from a Compound whose leaf particles are CG beads."""
    particles = list(compound.particles())
    if not particles:
        raise ValueError(
            "Compound has no particles; cannot build a CGsmiles graph."
        )
    index_of = {particle: index for index, particle in enumerate(particles)}
    meta_graph = nx.Graph()
    for index, particle in enumerate(particles):
        name = str(particle.name)
        meta_graph.add_node(
            index,
            fragname=_fragname(name, fragname_map),
            beadname=name,
            position=np.asarray(particle.pos, dtype=float),
        )
    for particle_a, particle_b in compound.bonds():
        meta_graph.add_edge(index_of[particle_a], index_of[particle_b], order=1)
    return meta_graph


def to_cgsmiles(system, fragname_map=None):
    """Write the coarse-grained level of a system as a CGsmiles string.

    The returned string describes the bead sequence and connectivity
    (including branches and rings) at the coarse-grained level.
    Append fragment definitions (e.g. ``"{#A=[>]CC[<]}"``)
    to obtain a fully resolvable CGsmiles string.

    Parameters
    ----------
    system : mbuild.path.Path or mbuild.Compound, required
        The coarse-grained system to convert. See ``to_cgsmiles_graph``.
    fragname_map : dict[str, str], optional
        Mapping of bead names to CGsmiles fragment names. Bead names not
        present in the map are used as fragment names directly.

    Returns
    -------
    str
        The CGsmiles graph string, e.g. ``"{[#A][#A]([#B][#B])[#A]}"``.
        Systems holding multiple disconnected molecules (e.g. a box of
        chains) are written as ``.``-separated segments, which CGsmiles
        reads as zero-order (non-bonded) connections.
    """
    from cgsmiles.write_cgsmiles import write_graph

    # TODO: cgsmiles' writer never emits the expansion operator (e.g.
    # "[#PEO]|10"); it is read-side only. Consider a post-processing pass
    # collapsing runs of consecutive identical [#name] tokens into [#name]|N.
    meta_graph = to_cgsmiles_graph(system, fragname_map=fragname_map)
    # cgsmiles' writer walks a single connected component, so write each
    # molecule separately and join with '.' (a zero-order connection)
    components = sorted(
        (sorted(component) for component in nx.connected_components(meta_graph)),
        key=lambda component: component[0],
    )
    segments = [
        write_graph(meta_graph.subgraph(component)) for component in components
    ]
    return "{" + ".".join(segments) + "}"
