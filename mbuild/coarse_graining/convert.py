"""Convert coarse-grained mBuild systems to CGsmiles graphs and strings.

A coarse-grained system is either an ``mbuild.path.Path`` or an
``mbuild.Compound`` whose leaf particles are the CG beads. Both carry
the same information (bead names, bead coordinates, and a bond graph)
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
    """Map a bead name to its CGsmiles fragment name via ``fragname_map``,
    returning the name unchanged when there is no mapping for it.
    """
    if fragname_map is not None:
        return fragname_map.get(name, name)
    return name


def _build_meta_graph(nodes, edges, fragname_map):
    """Assemble a CGsmiles meta graph from bead and bond iterables.

    Parameters
    ----------
    nodes : iterable of (int, str, array-like)
        One tuple per bead: its node index, bead name, and position.
    edges : iterable of (int, int)
        Pairs of bead indices that share a CG bond.
    fragname_map : dict[str, str] or None
        Mapping of bead names to CGsmiles fragment names.

    Returns
    -------
    networkx.Graph
        Each node carries ``fragname``, ``beadname``, and ``position``;
        each edge carries ``order=1``.
    """
    meta_graph = nx.Graph()
    for index, name, position in nodes:
        meta_graph.add_node(
            int(index),
            fragname=_fragname(name, fragname_map),
            beadname=name,
            position=np.asarray(position, dtype=float),
        )
    for u, v in edges:
        meta_graph.add_edge(int(u), int(v), order=1)
    return meta_graph


def _sorted_components(graph):
    """Connected components of a graph, each sorted and ordered by first node.

    Returns
    -------
    list of list of int
        Each inner list holds the sorted node ids of one connected
        component; components are ordered by their smallest node id.
    """
    return sorted(
        (sorted(component) for component in nx.connected_components(graph)),
        key=lambda component: component[0],
    )


def _path_meta_graph(path, fragname_map):
    """Meta graph from a Path's beads, coordinates, and bond graph."""
    # TODO: If we can set bond order in Path, grab it and set it here
    if len(path.bond_graph) != len(path.coordinates):
        path._extend_bond_graph()
    nodes = (
        (node, str(path.beads[node]), path.coordinates[node])
        for node in sorted(path.bond_graph.nodes)
    )
    return _build_meta_graph(nodes, path.bond_graph.edges, fragname_map)


def _compound_meta_graph(compound, fragname_map):
    """Meta graph from a Compound whose leaf particles are CG beads."""
    particles = list(compound.particles())
    if not particles:
        raise ValueError(
            "Compound has no particles; cannot build a CGsmiles graph."
        )
    index_of = {particle: index for index, particle in enumerate(particles)}
    nodes = (
        (index, str(particle.name), particle.pos)
        for index, particle in enumerate(particles)
    )
    edges = (
        (index_of[particle_a], index_of[particle_b])
        for particle_a, particle_b in compound.bonds()
    )
    return _build_meta_graph(nodes, edges, fragname_map)


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
    segments = [
        write_graph(meta_graph.subgraph(component))
        for component in _sorted_components(meta_graph)
    ]
    return "{" + ".".join(segments) + "}"
