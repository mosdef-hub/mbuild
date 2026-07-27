"""Backmap coarse-grained Paths to atomistic Compounds using CGsmiles.

This module bridges ``mbuild.path.Path`` objects (a bond graph with
coarse-grained coordinates and bead names) and the CGsmiles library
(https://github.com/gruenewald-lab/CGsmiles).

Example:

1. Build a CG path in mBuild.
2. Convert the path's bond graph to a CGsmiles-compatible meta graph
   where bead names become fragment names (``path_to_cgsmiles_graph``).
3. Resolve the meta graph to a molecular (atomistic) graph using
   CGsmiles fragment strings with SMILES + bonding descriptors
   (``resolve_path``).
4. Convert the resolved molecular graph back to an atomistic
   ``mbuild.Compound``, placing atoms according to the CG bead
   coordinates (``backmap_path``).

Unlike ``Polymer.build_from_path``, this approach works for any bond
graph topology, including beads with more than two edges (branch
points, stars, and crosslinks).

This utilizes the CGSmiles library (https://cgsmiles.readthedocs.io/en/stable/).
If you use this feature please cite the following paper:

CGsmiles: A Versatile Line Notation for Molecular Representations across Multiple Resolutions
Fabian Grünewald, Leif Seute, Riccardo Alessandri, Melanie König, and Peter C. Kroon
Journal of Chemical Information and Modeling 2025 65 (7), 3405-3419
DOI: 10.1021/acs.jcim.5c00064
"""

import logging
from itertools import combinations

import networkx as nx
import numpy as np

from mbuild.compound import Compound
from mbuild.utils.geometry import kabsch

logger = logging.getLogger(__name__)

__all__ = [
    "path_to_cgsmiles_graph",
    "path_to_cgsmiles",
    "compound_to_fragment_graph",
    "resolve_path",
    "backmap_path",
]


def path_to_cgsmiles_graph(path, fragname_map=None):
    """Convert a Path's bond graph to a CGsmiles-compatible meta graph.

    Each node in the returned graph carries the ``fragname`` attribute
    taken from the path's bead names, plus a ``position`` attribute taken from the path
    coordinates. All edges are given ``order=1``.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The coarse-grained path to convert.
    fragname_map : dict[str, str], optional
        Mapping of bead names to CGsmiles fragment names. Bead names not
        present in the map are used as fragment names directly. Useful
        when path bead names differ from the fragment
        names used in the CGsmiles fragment string.

    Returns
    -------
    networkx.Graph
        Meta graph usable with ``cgsmiles.MoleculeResolver.from_graph``.
    """
    #TODO: If we can set bond order in Path, grab it and set it here
    if len(path.bond_graph) != len(path.coordinates):
        path._extend_bond_graph()
    meta_graph = nx.Graph()
    for node in sorted(path.bond_graph.nodes):
        name = str(path.beads[node])
        if fragname_map is not None:
            name = fragname_map.get(name, name)
        meta_graph.add_node(
            int(node),
            fragname=name,
            position=np.asarray(path.coordinates[node], dtype=float),
        )
    for u, v in path.bond_graph.edges:
        meta_graph.add_edge(int(u), int(v), order=1)
    return meta_graph


def path_to_cgsmiles(path, fragname_map=None):
    """Write the coarse-grained level of a Path as a CGsmiles string.

    The returned string describes the bead sequence and connectivity
    (including branches and rings) at the coarse-grained level.
    Append fragment definitions (e.g. ``"{#A=[>]CC[<]}"``)
    to obtain a fully resolvable CGsmiles string.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The coarse-grained path to convert.
    fragname_map : dict[str, str], optional
        Mapping of bead names to CGsmiles fragment names. Bead names not
        present in the map are used as fragment names directly. Useful
        when path bead names differ from the fragment
        names used in the CGsmiles fragment string.

    Returns
    -------
    str
        The CGsmiles graph string, e.g. ``"{[#A][#A]([#B][#B])[#A]}"``.
        Paths holding multiple disconnected molecules (e.g. a box of
        chains) are written as ``.``-separated segments, which CGsmiles
        reads as zero-order (non-bonded) connections.
    """
    from cgsmiles.write_cgsmiles import write_graph

    # TODO: cgsmiles' writer never emits the expansion operator (e.g.
    # "[#PEO]|10"); it is read-side only. Consider a post-processing pass
    # collapsing runs of consecutive identical [#name] tokens into [#name]|N.
    meta_graph = path_to_cgsmiles_graph(path, fragname_map=fragname_map)
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


_DESCRIPTOR_PREFIXES = ("<", ">", "$", "!")




def compound_to_fragment_graph(compound, fragname):
    """Convert a tagged mBuild Compound into a CGsmiles fragment graph.

    This lets an mBuild Compound *define* a fragment's atomistic
    structure, in place of a SMILES string in the CGsmiles fragment
    block. Bonding descriptors are read from particle tags: tag the
    connection atoms with descriptor strings, e.g.
    ``mb.load("C{$}C{$}", smiles=True)`` is equivalent to the fragment
    SMILES ``[$]CC[$]``, and ``"C{$,$}C{$}"`` puts two descriptors on
    the first carbon (a comma separates multiple descriptors on one
    atom). Descriptors follow the usual matching rules (``<`` bonds
    ``>``, ``$`` bonds ``$``, labels like ``$br`` must match).

    Hydrogen particles are folded into the heavy atoms' ``hcount``;
    the resolver removes hydrogens consumed by inter-fragment bonds
    automatically, so the compound should be saturated.

    Parameters
    ----------
    compound : mbuild.Compound, required
        The atomistic fragment. Connection atoms must carry descriptor
        tags (set via tagged SMILES or ``particle.particle_tag = "$"``).
    fragname : str, required
        The fragment name (must match the bead/fragment name it defines).

    Returns
    -------
    networkx.Graph
        A fragment graph with the same node/edge attributes as those
        returned by ``cgsmiles.read_fragments``, usable in a
        ``MoleculeResolver`` fragment dict.
    """
    heavy = [
        particle
        for particle in compound.particles()
        if _particle_symbol(particle) != "H"
    ]
    if not heavy:
        raise ValueError(
            f"Compound for fragment '{fragname}' has no heavy atoms."
        )
    index_of = {particle: index for index, particle in enumerate(heavy)}
    hcounts = dict.fromkeys(heavy, 0)

    heavy_bonds = []
    for particle_a, particle_b, data in compound.bonds(return_bond_order=True):
        symbol_a = _particle_symbol(particle_a)
        symbol_b = _particle_symbol(particle_b)
        if symbol_a == "H" or symbol_b == "H":
            if symbol_a != "H":
                hcounts[particle_a] += 1
            elif symbol_b != "H":
                hcounts[particle_b] += 1
            continue
        order = data.get("bond_order", data.get("bo"))
        if order is None or order == "default":
            order = 1
        heavy_bonds.append((particle_a, particle_b, float(order)))

    graph = nx.Graph()
    for particle in heavy:
        index = index_of[particle]
        attributes = {
            "element": _particle_symbol(particle),
            "charge": int(particle.charge or 0),
            "aromatic": False,
            "fragname": fragname,
            "fragid": 0,
            "weight": 1,
            "atomname": f"{_particle_symbol(particle)}{index}",
            "hcount": hcounts[particle],
        }
        descriptors = _descriptors_from_tag(particle.particle_tag)
        if descriptors:
            attributes["bonding"] = descriptors
        graph.add_node(index, **attributes)

    for particle_a, particle_b, order in heavy_bonds:
        if order == 1.5:
            graph.nodes[index_of[particle_a]]["aromatic"] = True
            graph.nodes[index_of[particle_b]]["aromatic"] = True
            order = 1.5
        else:
            order = int(order)
        graph.add_edge(index_of[particle_a], index_of[particle_b], order=order)
    return graph


def resolve_path(
    path, fragments=None, fragname_map=None, legacy=True, validate=True, templates=None
):
    """Resolve a coarse-grained Path to a molecular graph using CGsmiles.

    The path's bond graph is handed to CGsmiles as the lowest resolution
    graph and resolved through each fragment string down to the highest
    resolution (by default, all-atom).

    Parameters
    ----------
    path : mbuild.path.Path, required
        The coarse-grained path to resolve.
    fragments : str or list of str, optional
        CGsmiles fragment string(s) defining each resolution, e.g.
        ``"{#A=[>]CC[<],#B=[>]CC(C)[<]}"``. Fragment names must match
        the path's bead names (after applying ``fragname_map``), and the
        bonding descriptors (``[<]``, ``[>]``, ``[$]``) must provide
        enough connection points for the maximum degree of the beads.
        Multiple strings resolve through intermediate resolutions in
        order, with the last string interpreted as all-atom SMILES.
        May be omitted when ``templates`` defines every fragment.
    fragname_map : dict[str, str], optional
        Mapping of bead names to CGsmiles fragment names.
        See ``path_to_cgsmiles_graph``.
    templates : dict[str, mbuild.Compound], optional
        Tagged compounds that *define* fragments not present in the
        ``fragments`` string (see ``compound_to_fragment_graph``).
        Fragment names already defined in the string keep the string
        definition. With ``fragments=None``, the templates must define
        every fragment used by the path.
    legacy : bool, default True
        Bonding descriptor matching convention passed through to
        ``cgsmiles.MoleculeResolver`` (legacy=BigSmiles convention).
    validate : bool, default True
        If ``True``, verify that every bond in the path's bond graph was
        realized as an atomistic bond between the corresponding
        fragments. CGsmiles silently skips CG bonds when a fragment runs
        out of bonding descriptors, which would otherwise produce a
        disconnected molecule with hydrogen-capped stumps.

    Returns
    -------
    meta_graph : networkx.Graph
        The second to last resolution graph with each node annotated
        with its fragment ``graph``.
    molecule : networkx.Graph
        The fully resolved (atomistic) molecule graph. Nodes carry
        pysmiles-style attributes (``element``, ``charge``, ...) and
        edges carry bond ``order``.
    node_to_beads : dict[int, list[int]]
        Maps each atom (node in ``molecule``) to the sorted list of
        path bead indices it descends from. Atoms claimed by multiple
        beads (via the squash operator) map to multiple indices.
    """
    import re

    from cgsmiles.resolve import MoleculeResolver

    if fragments is None and not templates:
        raise ValueError(
            "Provide a CGsmiles fragment string, fragment-defining "
            "templates, or both."
        )
    if isinstance(fragments, str):
        fragments = [fragments]

    if fragments:
        elements = re.findall(r"\{[^\}]+\}", ".".join(fragments))
        fragment_dicts = MoleculeResolver.read_fragment_strings(
            elements, last_all_atom=True
        )
    else:
        fragment_dicts = [{}]

    # compounds define any fragments the strings left undefined at the
    # final all-atom resolution
    if templates:
        final_dict = fragment_dicts[-1]
        for fragname, compound in templates.items():
            if fragname not in final_dict:
                final_dict[fragname] = compound_to_fragment_graph(compound, fragname)

    meta_graph = path_to_cgsmiles_graph(path, fragname_map=fragname_map)
    resolver = MoleculeResolver(
        molecule_graph=meta_graph,
        fragment_dicts=fragment_dicts,
        last_all_atom=True,
        legacy=legacy,
    )

    # Track which original path bead(s) each node descends from across
    # resolutions; each resolution's `fragid` refers to nodes of the
    # previous resolution graph.
    node_to_beads = None
    for meta_graph, molecule in resolver.resolve_iter():
        new_map = {}
        for node, data in molecule.nodes(data=True):
            fragids = data.get("fragid", [])
            if node_to_beads is None:
                beads = set(fragids)
            else:
                beads = set()
                for fragid in fragids:
                    beads.update(node_to_beads[fragid])
            new_map[node] = sorted(beads)
        node_to_beads = new_map

    if validate:
        _validate_connectivity(path, molecule, node_to_beads)

    return meta_graph, molecule, node_to_beads




def backmap_path(
    path,
    fragments=None,
    fragname_map=None,
    templates=None,
    placement="auto",
    jitter=0.05,
    energy_minimize=False,
    seed=42,
    legacy=True,
    validate=True,
):
    """Backmap a coarse-grained Path to an atomistic mbuild Compound.

    Resolves the path's bond graph to an atomistic molecule via CGsmiles
    (see ``resolve_path``) and builds a Compound whose hierarchy mirrors
    the path: one child Compound per CG bead, named after the bead,
    containing the atoms that bead resolved to. Atoms are positioned
    based on the CG bead coordinates so the compound retains the path's
    conformation.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The coarse-grained path to backmap.
    fragments : str or list of str, optional
        CGsmiles fragment string(s) defining each resolution, e.g.
        ``"{#A=[>]CC[<],#B=[>]CC(C)[<]}"``. Fragment names must match
        the path's bead names (after applying ``fragname_map``), and the
        bonding descriptors (``[<]``, ``[>]``, ``[$]``) must provide
        enough connection points for the maximum degree of the beads.
        Multiple strings resolve through intermediate resolutions in
        order, with the last string interpreted as all-atom SMILES.
        May be omitted when ``templates`` defines every fragment.
    fragname_map : dict[str, str], optional
        Mapping of bead names to CGsmiles fragment names.
        See ``path_to_cgsmiles_graph``.
    templates : dict[str, mbuild.Compound], optional
        Atomistic template compounds keyed by fragment name. Templates
        play two roles:

        1. **Structure**: a template whose fragment is *not* defined in
           the ``fragments`` string defines the fragment itself — its
           atoms, bonds, and (via particle tags carrying bonding
           descriptors, e.g. ``mb.load("C{$}C{$}", smiles=True)``) its
           connection points. No SMILES needed; see
           ``compound_to_fragment_graph``. With ``fragments=None`` the
           templates must define every fragment.
        2. **Geometry**: every templated bead takes its local
           coordinates (and hydrogen positions) from the template
           instead of an RDKit embedding. For fragments defined by a
           SMILES string, template heavy atoms are matched by fragment
           atom order (element-checked) or explicitly, by tagging every
           heavy atom with its fragment atom index, e.g.
           ``mb.load("O{1}(C{0})C{2}", smiles=True)`` for the fragment
           SMILES ``COC``. Fragments defined by the compound itself
           match by construction.

        The template must be at least as saturated as the resolved
        fragment (extra hydrogens, e.g. at junctions, are ignored).
        Beads without a template fall back to the ``placement``
        behavior.
    placement : str, default "auto"
        How to generate initial atom positions for beads not covered by
        ``templates``.

        - ``"rdkit"``: embed each bead's fragment with RDKit to obtain
          valid local geometry (identical fragments are embedded once
          and reused), rotate it so its junction atoms point toward the
          neighboring beads, and center it on the bead coordinate. Cost
          scales with the number of unique fragments, not molecule size.
        - ``"jitter"``: place each atom at its bead coordinate plus a
          small random displacement of scale ``jitter``.
        - ``"auto"``: use ``"rdkit"`` if RDKit is available and the
          embedding succeeds, otherwise fall back to ``"jitter"``.
        - ``"template"``: require a template for every fragment and
          raise if one is missing.

        Template-placed fragments are rotated toward their neighboring
        beads exactly like RDKit-placed ones.
    jitter : float, default 0.05 (nm)
        Scale of the random displacement used by ``"jitter"`` placement.
    energy_minimize : bool, default False
        If ``True``, relax the backmapped compound with
        ``mbuild.simulation.energy_minimize`` (OpenBabel UFF). For large
        molecules the hoomd-based methods in ``mbuild.simulation`` are
        significantly faster; leave this ``False`` and pass the returned
        compound to one of those instead.
    seed : int, default 42
        Random seed used for jitter displacements and RDKit embedding.
    legacy : bool, default True
        Bonding descriptor matching convention. See ``resolve_path``.
    validate : bool, default True
        Verify that every CG bond was realized as an atomistic bond.
        See ``resolve_path``.

    Returns
    -------
    mbuild.Compound
        The atomistic compound. Bonds between atoms follow the resolved
        molecular graph, including the bonds between fragments dictated
        by the path's bond graph. If the path holds several disconnected
        molecules (e.g. a box of chains), the hierarchy gains an
        intermediate level: root -> one ``Molecule`` compound per
        connected component -> bead compounds; a single-molecule path
        keeps its bead compounds directly under the root.

    Notes
    -----
    Initial configurations produced by this method have distorted bonds
    at the fragment-fragment connections (and overlapping atoms when
    using ``"jitter"`` placement); an energy minimization is required
    before using them in a simulation.

    See mbuild.Simulation for minimization methods.
    """
    meta_graph, molecule, node_to_beads = resolve_path(
        path,
        fragments,
        fragname_map=fragname_map,
        legacy=legacy,
        validate=validate,
        templates=templates,
    )

    atom_nodes = list(molecule.nodes)
    bead_positions = {
        int(node): np.asarray(path.coordinates[node], dtype=float)
        for node in path.bond_graph.nodes
    }

    # Anchor of each atom: mean position of the bead it descends from
    anchors = {}
    for node in atom_nodes:
        beads = node_to_beads[node]
        if not beads:
            raise ValueError(
                f"Atom {node} does not trace back to any path bead; "
                "cannot assign a position."
            )
        anchors[node] = np.mean([bead_positions[b] for b in beads], axis=0)

    if placement not in ("auto", "rdkit", "jitter", "template"):
        raise ValueError(
            f"Argument {placement=} is invalid. "
            "Pass one of 'auto', 'rdkit', 'jitter', 'template'."
        )
    if placement == "template" and not templates:
        raise ValueError(
            "placement='template' requires the `templates` argument."
        )

    bead_fragnames = {}
    for node in path.bond_graph.nodes:
        name = str(path.beads[node])
        if fragname_map is not None:
            name = fragname_map.get(name, name)
        bead_fragnames[int(node)] = name

    positions = _generate_positions(
        molecule=molecule,
        node_to_beads=node_to_beads,
        anchors=anchors,
        bead_fragnames=bead_fragnames,
        placement=placement,
        templates=templates,
        jitter=jitter,
        seed=seed,
    )

    compound = _molecule_to_compound(
        path, molecule, node_to_beads, positions, fragname_map
    )

    if energy_minimize:
        from mbuild.simulation import energy_minimize as e_min

        e_min(compound)
    return compound


def _molecule_to_compound(path, molecule, node_to_beads, positions, fragname_map):
    """Build an mbuild Compound from a resolved molecule graph.

    The hierarchy is one child Compound per CG bead (named after the
    bead) containing that bead's atoms. Atoms descending from several
    beads are placed in the compound of the lowest bead index. If the
    path's bond graph holds multiple disconnected molecules, an
    intermediate ``Molecule`` compound groups each molecule's beads.
    """
    from ele import element_from_symbol
    from ele.exceptions import ElementError

    root = Compound(name="Backmapped")
    bead_compounds = {}
    for node in sorted(path.bond_graph.nodes):
        bead_compounds[int(node)] = Compound(name=str(path.beads[node]))

    particles = {}
    for node, data in molecule.nodes(data=True):
        symbol = data.get("element")
        try:
            element = element_from_symbol(symbol)
        except (ElementError, TypeError):
            element = None
        particle = Compound(
            name=symbol,
            pos=positions[node],
            element=element,
            charge=data.get("charge", 0) or None,
        )
        particles[node] = particle
        bead_compounds[node_to_beads[node][0]].add(particle)

    components = sorted(
        (
            sorted(int(node) for node in component)
            for component in nx.connected_components(path.bond_graph)
        ),
        key=lambda component: component[0],
    )
    if len(components) == 1:
        for bead_index in sorted(bead_compounds):
            bead_compound = bead_compounds[bead_index]
            if bead_compound.children:
                root.add(bead_compound)
    else:
        for component in components:
            molecule_compound = Compound(name="Molecule")
            for bead_index in component:
                bead_compound = bead_compounds[bead_index]
                if bead_compound.children:
                    molecule_compound.add(bead_compound)
            root.add(molecule_compound, label="Molecule[$]")

    for u, v, data in molecule.edges(data=True):
        root.add_bond((particles[u], particles[v]), bond_order=float(data["order"]))
    return root


def _generate_positions(
    molecule,
    node_to_beads,
    anchors,
    bead_fragnames,
    placement,
    templates,
    jitter,
    seed,
):
    """Compute a position for every atom in the resolved molecule.

    For each bead, local fragment geometry is taken from a user template
    (if one covers the bead's fragment), or an RDKit embedding of the
    fragment (cached per unique fragment), or jitter around the bead
    anchor. Template- and RDKit-placed fragments are rotated so the
    atoms bonded to other beads point toward those beads, then centered
    on the bead anchor.
    """
    rdkit_available = True
    if placement in ("auto", "rdkit"):
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            if placement == "rdkit":
                raise ImportError(
                    "placement='rdkit' requires the `rdkit` package. "
                    "Install it with `conda install -c conda-forge rdkit`."
                )
            rdkit_available = False
            logger.warning(
                "RDKit not available; falling back to jitter placement."
            )

    # Group atoms by primary bead
    bead_to_atoms = {}
    for node in molecule.nodes:
        bead_to_atoms.setdefault(node_to_beads[node][0], []).append(node)

    rng = np.random.default_rng(seed)
    embed_cache = {}
    positions = {}
    for bead, atoms in sorted(bead_to_atoms.items()):
        atoms = sorted(atoms)
        atom_set = set(atoms)
        rel_index = {node: i for i, node in enumerate(atoms)}
        fragname = bead_fragnames.get(bead)

        local_xyz = None
        if templates is not None and fragname in templates:
            local_xyz = _template_local_coords(
                molecule, atoms, templates[fragname], fragname
            )
        elif placement == "template":
            raise ValueError(
                f"placement='template' but no template was provided for "
                f"fragment '{fragname}'."
            )
        elif placement in ("auto", "rdkit") and rdkit_available:
            try:
                local_xyz = _embedded_fragment_coords(
                    molecule, atoms, rel_index, embed_cache, seed
                )
            except Exception as error:
                if placement == "rdkit":
                    raise
                logger.warning(
                    "RDKit embedding failed for fragment '%s' (%s); "
                    "falling back to jitter placement for this bead.",
                    fragname,
                    error,
                )

        if local_xyz is None:  # jitter placement or rdkit fallback
            for node in atoms:
                positions[node] = anchors[node] + rng.normal(scale=jitter, size=3)
            continue

        group_anchor = np.mean([anchors[node] for node in atoms], axis=0)

        # Rotate the fragment so its junction atoms (those bonded to
        # atoms in other beads) point toward the neighboring beads
        sources = []
        targets = []
        for node in atoms:
            for neighbor in molecule.neighbors(node):
                if neighbor in atom_set:
                    continue
                source = local_xyz[rel_index[node]]
                target = anchors[neighbor] - group_anchor
                source_norm = np.linalg.norm(source)
                target_norm = np.linalg.norm(target)
                if source_norm > 1e-8 and target_norm > 1e-8:
                    sources.append(source / source_norm)
                    targets.append(target / target_norm)
        if sources:
            rotation = kabsch(np.array(sources), np.array(targets))
            local_xyz = local_xyz @ rotation.T

        for node in atoms:
            positions[node] = group_anchor + local_xyz[rel_index[node]]
    return positions


def _embedded_fragment_coords(molecule, atoms, rel_index, embed_cache, seed):
    """Local fragment coordinates from a cached RDKit embedding.

    Fragment instances created from the same CGsmiles fragment have
    their atoms in the same relative node order, so element, charge and
    relative-edge tuples identify structurally identical fragments and
    index-based reuse of the embedded coordinates is safe.
    """
    elements = tuple(molecule.nodes[n].get("element", "*") for n in atoms)
    charges = tuple(molecule.nodes[n].get("charge", 0) for n in atoms)
    edges = tuple(
        sorted(
            (
                min(rel_index[u], rel_index[v]),
                max(rel_index[u], rel_index[v]),
                data.get("order", 1),
            )
            for u, v, data in molecule.subgraph(atoms).edges(data=True)
        )
    )
    key = (elements, charges, edges)
    if key not in embed_cache:
        embed_cache[key] = _embed_fragment(elements, charges, edges, seed)
    return embed_cache[key]


def _tag_as_index(particle):
    """Interpret a particle's tag as a fragment atom index, if possible.

    Non-integer tags (e.g. the "head"/"tail" tags used by the Polymer
    workflow) are ignored so templates carrying them still match by
    atom order.
    """
    tag = particle.particle_tag
    if tag is None:
        return None
    try:
        return int(str(tag))
    except ValueError:
        return None


def _template_local_coords(molecule, atoms, template, fragname):
    """Local coordinates for a bead's atoms taken from an mBuild template.

    Heavy atoms in the resolved molecule are matched to the template's
    heavy particles either explicitly — every heavy particle tagged with
    its fragment atom index, e.g. ``mb.load("O{1}(C{0})C{2}",
    smiles=True)`` — or, for untagged templates, by fragment atom order
    (the CGsmiles ``mapping`` node attribute), in which case the
    template's heavy atoms must appear in the same order as in the
    fragment SMILES. Each hydrogen is assigned the position of an unused
    hydrogen bonded to its parent's matched template particle; surplus
    template hydrogens (e.g. saturating a junction) are ignored. Returns
    coordinates in nm, centered on the centroid of the used atoms.
    """

    def symbol(particle):
        if particle.element is not None:
            return particle.element.symbol
        return particle.name

    heavy_nodes = []
    for node in atoms:
        if molecule.nodes[node].get("element") == "H":
            continue
        mapping = molecule.nodes[node].get("mapping")
        if not mapping or len(mapping) != 1:
            raise ValueError(
                f"Cannot use a template for fragment '{fragname}': atom "
                f"{node} has no unique fragment mapping (squashed fragments "
                "are not supported by template placement)."
            )
        heavy_nodes.append(node)
    heavy_nodes.sort(key=lambda n: molecule.nodes[n]["mapping"][0][1])

    template_heavy = [p for p in template.particles() if symbol(p) != "H"]
    if len(template_heavy) != len(heavy_nodes):
        raise ValueError(
            f"Template for fragment '{fragname}' has {len(template_heavy)} "
            f"heavy atoms but the resolved fragment has {len(heavy_nodes)}."
        )

    # Explicit correspondence via particle tags: heavy atoms tagged with
    # their fragment atom index (e.g. mb.load("O{1}(C{0})C{2}")) are
    # matched by tag; untagged templates are matched by atom order.
    tag_indices = [_tag_as_index(particle) for particle in template_heavy]
    n_tagged = sum(1 for tag in tag_indices if tag is not None)
    if 0 < n_tagged < len(template_heavy):
        raise ValueError(
            f"Template for fragment '{fragname}' has fragment-index tags on "
            f"{n_tagged} of {len(template_heavy)} heavy atoms; tag all heavy "
            "atoms or none."
        )
    if n_tagged:
        if len(set(tag_indices)) != len(tag_indices):
            raise ValueError(
                f"Template for fragment '{fragname}' has duplicate "
                "fragment-index tags on its heavy atoms."
            )
        particle_by_index = dict(zip(tag_indices, template_heavy))
        pairs = []
        for node in heavy_nodes:
            local_index = molecule.nodes[node]["mapping"][0][1]
            particle = particle_by_index.get(local_index)
            if particle is None:
                raise ValueError(
                    f"Template for fragment '{fragname}' has no heavy atom "
                    f"tagged {{{local_index}}}; the fragment uses atom "
                    f"indices {sorted(molecule.nodes[n]['mapping'][0][1] for n in heavy_nodes)}."
                )
            pairs.append((node, particle))
    else:
        pairs = list(zip(heavy_nodes, template_heavy))

    coords = {}
    matched = {}
    for position, (node, particle) in enumerate(pairs):
        element = molecule.nodes[node]["element"]
        if element != symbol(particle):
            hint = (
                "Check the fragment-index tags."
                if n_tagged
                else "Template heavy atoms must appear in the same order as "
                "in the fragment SMILES, or be tagged with fragment atom "
                "indices (e.g. mb.load(\"O{1}(C{0})C{2}\", smiles=True))."
            )
            raise ValueError(
                f"Template for fragment '{fragname}' does not match: heavy "
                f"atom {position} is '{symbol(particle)}' but the fragment "
                f"expects '{element}'. {hint}"
            )
        coords[node] = np.asarray(particle.pos, dtype=float)
        matched[node] = particle

    # hydrogens available on each matched template heavy atom
    available_hydrogens = {particle: [] for particle in template_heavy}
    for particle_a, particle_b in template.bonds():
        if symbol(particle_a) == "H" and particle_b in available_hydrogens:
            available_hydrogens[particle_b].append(particle_a)
        elif symbol(particle_b) == "H" and particle_a in available_hydrogens:
            available_hydrogens[particle_a].append(particle_b)

    atom_set = set(atoms)
    for node in atoms:
        if molecule.nodes[node].get("element") != "H":
            continue
        parents = [
            neighbor
            for neighbor in molecule.neighbors(node)
            if neighbor in atom_set and neighbor in matched
        ]
        if not parents:
            raise ValueError(
                f"Cannot use a template for fragment '{fragname}': hydrogen "
                f"atom {node} is not bonded to a heavy atom in the fragment."
            )
        pool = available_hydrogens[matched[parents[0]]]
        if not pool:
            raise ValueError(
                f"Template for fragment '{fragname}' has too few hydrogens "
                f"on heavy atom {heavy_nodes.index(parents[0])}; the template "
                "must be at least as saturated as the resolved fragment."
            )
        coords[node] = np.asarray(pool.pop(0).pos, dtype=float)

    local_xyz = np.array([coords[node] for node in atoms])
    return local_xyz - local_xyz.mean(axis=0)


def _embed_fragment(elements, charges, edges, seed):
    """Embed a single fragment with RDKit.

    The fragment is described by per-atom element symbols and formal
    charges plus (index, index, order) edge tuples. Atoms bonded to
    other fragments keep an open valence; RDKit treats them as radicals,
    which does not disturb the embedded geometry. Returns coordinates in
    nm, centered on the fragment centroid.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    bond_type_map = {
        0: Chem.BondType.ZERO,
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        4: Chem.BondType.QUADRUPLE,
        1.5: Chem.BondType.AROMATIC,
    }

    fragment = Chem.RWMol()
    for element, charge in zip(elements, charges):
        atom = Chem.Atom(element)
        atom.SetFormalCharge(charge)
        atom.SetNoImplicit(True)  # all H atoms are explicit after resolving
        fragment.AddAtom(atom)
    for u, v, order in edges:
        fragment.AddBond(u, v, bond_type_map.get(order, Chem.BondType.SINGLE))
        if order == 1.5:
            fragment.GetAtomWithIdx(u).SetIsAromatic(True)
            fragment.GetAtomWithIdx(v).SetIsAromatic(True)

    fragment = fragment.GetMol()
    Chem.SanitizeMol(fragment)
    if AllChem.EmbedMolecule(fragment, randomSeed=seed) != 0:
        if AllChem.EmbedMolecule(fragment, randomSeed=seed, useRandomCoords=True) != 0:
            raise ValueError(f"RDKit failed to embed fragment {elements}.")
    try:
        AllChem.UFFOptimizeMolecule(fragment)
    except Exception:  # UFF parameters can be missing; keep raw DG coordinates
        pass

    # Angstrom -> nm
    xyz = fragment.GetConformer().GetPositions() / 10.0
    return xyz - xyz.mean(axis=0)


def _validate_connectivity(path, molecule, node_to_beads):
    """Check that every CG bond was realized as an atomistic bond.

    CGsmiles skips a CG bond without raising when it cannot find a pair
    of matching bonding descriptors between the two fragments, leaving a
    disconnected molecule behind. Raises a ValueError naming the
    unrealized CG bonds.
    """
    realized = set()
    for u, v in molecule.edges:
        for bead_u in node_to_beads[u]:
            for bead_v in node_to_beads[v]:
                if bead_u != bead_v:
                    realized.add((min(bead_u, bead_v), max(bead_u, bead_v)))
    # atoms shared by several beads (squash operator) also realize bonds
    for beads in node_to_beads.values():
        for bead_u, bead_v in combinations(beads, 2):
            realized.add((min(bead_u, bead_v), max(bead_u, bead_v)))

    expected = {
        (min(int(u), int(v)), max(int(u), int(v))) for u, v in path.bond_graph.edges
    }
    missing = sorted(expected - realized)
    if missing:
        details = ", ".join(
            f"({u} '{path.beads[u]}' — {v} '{path.beads[v]}')" for u, v in missing[:10]
        )
        if len(missing) > 10:
            details += ", ..."
        raise ValueError(
            f"{len(missing)} CG bond(s) in the path were not realized as "
            f"atomistic bonds: {details}. This usually means a fragment ran "
            "out of bonding descriptors; each fragment needs at least as many "
            "descriptors as the highest degree of the beads that use it. "
            "Pass validate=False to skip this check."
        )


def _particle_symbol(particle):
    """Element symbol of a particle, falling back to its name."""
    if particle.element is not None:
        return particle.element.symbol
    return particle.name


def _descriptors_from_tag(tag):
    """Parse CGsmiles bonding descriptors out of a particle tag.

    Tags are comma-separated; tokens starting with a descriptor prefix
    (``<``, ``>``, ``$``, ``!``) are returned with the bond order suffix
    CGsmiles expects (``"$br"`` -> ``"$br1"``). Other tokens (e.g.
    ``head``/``tail`` or integer correspondence tags) are ignored.
    """
    if tag is None:
        return []
    descriptors = []
    for token in str(tag).split(","):
        token = token.strip()
        if token and token[0] in _DESCRIPTOR_PREFIXES:
            descriptors.append(token if token[-1].isdigit() else token + "1")
    return descriptors
