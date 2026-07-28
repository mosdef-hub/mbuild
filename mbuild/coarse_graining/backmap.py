"""Backmap coarse-grained systems to atomistic Compounds using CGsmiles.

The CG system can be an ``mbuild.path.Path`` or an ``mbuild.Compound``
whose leaf particles are the beads. It is handed to CGsmiles as the
lowest-resolution graph and resolved to molecular detail:

This process utilizes the CGSmiles python package.
If you use this feature please cite:

CGsmiles: A Versatile Line Notation for Molecular Representations
across Multiple Resolutions. Fabian Grünewald, Leif Seute, Riccardo
Alessandri, Melanie König, and Peter C. Kroon. Journal of Chemical
Information and Modeling 2025 65 (7), 3405-3419.
DOI: 10.1021/acs.jcim.5c00064
"""

import logging
from itertools import combinations

import networkx as nx
import numpy as np

from mbuild.compound import Compound

from .convert import _sorted_components, to_cgsmiles_graph
from .fragments import compound_to_fragment_graph
from .placement import generate_positions

logger = logging.getLogger(__name__)

__all__ = ["resolve", "backmap"]


def resolve(
    system, fragments=None, fragname_map=None, legacy=True, validate=True, templates=None
):
    """Resolve a coarse-grained system to a molecular graph using CGsmiles.

    The system's bond graph is handed to CGsmiles as the lowest
    resolution graph and resolved through each fragment string down to
    the highest resolution (by default, all-atom).

    Parameters
    ----------
    system : mbuild.path.Path or mbuild.Compound, required
        The coarse-grained system to resolve. For a ``Compound``, each
        leaf particle is treated as one CG bead (see
        ``to_cgsmiles_graph``).
    fragments : str or list of str, optional
        CGsmiles fragment string(s) defining each resolution, e.g.
        ``"{#A=[>]CC[<],#B=[>]CC(C)[<]}"``. Fragment names must match
        the bead names (after applying ``fragname_map``), and the
        bonding descriptors (``[<]``, ``[>]``, ``[$]``) must provide
        enough connection points for the maximum degree of the beads.
        Multiple strings resolve through intermediate resolutions in
        order, with the last string interpreted as all-atom SMILES.
        May be omitted when ``templates`` defines every fragment.
    fragname_map : dict[str, str], optional
        Mapping of bead names to CGsmiles fragment names.
        See ``to_cgsmiles_graph``.
    templates : dict[str, mbuild.Compound], optional
        Template mBuild Compounds keyed by fragment name. A template
        does two things for the beads that use it:

        1. Defines the fragment's structure, when its name is not in
           the ``fragments`` string. The compound's particles and bonds
           become the fragment's atoms and bonds, and particle tags
           mark the connection atoms: ``mb.load("C{$}C{$}",
           smiles=True)`` is equivalent to the fragment SMILES
           ``[$]CC[$]``, or set ``particle.particle_tag = "$"`` on a
           hand-built compound.
        2. Supplies the bead's local geometry. The template's
           coordinates (including hydrogen positions) are used in
           place of an RDKit embedding.
    legacy : bool, default True
        Bonding descriptor matching convention passed through to
        ``cgsmiles.MoleculeResolver`` (legacy=BigSmiles convention).
    validate : bool, default True
        If ``True``, verify that every bond in the CG bond graph was
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
        CG bead indices it descends from. Atoms claimed by multiple
        beads (via the squash operator) map to multiple indices.
    """
    cg_graph = to_cgsmiles_graph(system, fragname_map=fragname_map)
    return _resolve_graph(
        cg_graph,
        fragments,
        legacy=legacy,
        validate=validate,
        templates=templates,
    )


def _resolve_graph(cg_graph, fragments, legacy, validate, templates):
    """Resolve a prebuilt CG meta graph. See ``resolve``."""
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

    resolver = MoleculeResolver(
        molecule_graph=cg_graph,
        fragment_dicts=fragment_dicts,
        last_all_atom=True,
        legacy=legacy,
    )

    # Track which original CG bead(s) each node descends from across
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
        _validate_connectivity(cg_graph, molecule, node_to_beads)

    return meta_graph, molecule, node_to_beads


def backmap(
    system,
    fragments=None,
    fragname_map=None,
    templates=None,
    placement="rdkit",
    seed=42,
    legacy=True,
    validate=True,
):
    """Backmap a coarse-grained system to an atomistic mbuild Compound.

    Resolves the system's bond graph to a higher resolution structure via
    CGsmiles (see ``resolve``) and builds a Compound whose hierarchy
    mirrors the CG structure: one child Compound per CG bead, named
    after the bead, containing the atoms that bead resolved to. Atoms
    are positioned based on the CG bead coordinates so the compound
    retains the CG conformation.

    Parameters
    ----------
    system : mbuild.path.Path or mbuild.Compound, required
        The coarse-grained system to backmap. For a ``Compound``, each
        leaf particle is treated as one CG bead: particle names become
        fragment names, particle positions become bead positions, and
        the compound's bonds become the CG bond graph (see
        ``to_cgsmiles_graph``).
    fragments : str or list of str, optional
        CGsmiles fragment string(s) defining each resolution, e.g.
        ``"{#A=[>]CC[<],#B=[>]CC(C)[<]}"``. Fragment names must match
        the bead names in the Compound or Path, and the
        bonding descriptors (``[<]``, ``[>]``, ``[$]``) must provide
        enough connection points for the maximum degree of the CG beads.
        Multiple strings resolve through intermediate resolutions in
        order, with the last string interpreted as the highest resolution.
        May be omitted when ``templates`` defines every fragment.
    fragname_map : dict[str, str], optional
        Mapping of bead names to CGsmiles fragment names.
        See ``to_cgsmiles_graph``.
    templates : dict[str, mbuild.Compound], optional
        Template mBuild Compounds keyed by fragment name. A template
        does two things for the beads that use it:

        1. Defines the fragment's structure, when its name is not in
           the ``fragments`` string. The compound's particles and bonds
           become the fragment's atoms and bonds, and particle tags
           mark the connection atoms: ``mb.load("C{$}C{$}",
           smiles=True)`` is equivalent to the fragment SMILES
           ``[$]CC[$]``, or set ``particle.particle_tag = "$"`` on a
           hand-built compound.
        2. Supplies the bead's local geometry. The template's
           coordinates (including hydrogen positions) are used in
           place of an RDKit embedding.
    placement : str, default "rdkit"
        How to generate initial atom positions for beads not covered by
        ``templates``.

        - ``"rdkit"``: embed each bead's fragment with RDKit to obtain
          valid local geometry (identical fragments are embedded once
          and reused), rotated it so its junction atoms point toward the
          neighboring beads, and center it on the bead coordinate. Cost
          scales with the number of unique fragments, not molecule
          size. Raises if an embedding fails; provide a template for
          that fragment instead.
        - ``"template"``: require a template for every fragment and
          raise if one is missing. No RDKit needed.

        Template-placed fragments are rotated toward their neighboring
        beads exactly like RDKit-placed ones.
    seed : int, default 42
        Random seed used for RDKit embedding.
    legacy : bool, default True
        Bonding descriptor matching convention. See ``resolve``.
    validate : bool, default True
        Verify that every CG bond was realized as an atomistic bond.
        See ``resolve``.

    Returns
    -------
    mbuild.Compound
        The atomistic compound. Bonds between atoms follow the resolved
        molecular graph, including the bonds between fragments dictated
        by the CG bond graph. If the system holds several disconnected
        molecules (e.g. a box of chains), the hierarchy gains an
        intermediate level: root -> one ``Molecule`` compound per
        connected component -> bead compounds; a single-molecule system
        keeps its bead compounds directly under the root.

    Notes
    -----
    Initial configurations produced by this method have distorted bonds
    at the fragment-fragment connections; an energy minimization is
    required before using them in a simulation.

    See mbuild.Simulation for minimization methods.

    Examples
    --------
    1. CGsmiles string only: The string defines everything and RDKit
       generates the geometry::

           path = straight_line(spacing=0.35, N=8, bead_name="PEO")
           compound = path.backmap("{#PEO=[>]COC[<]}")

    2. mBuild compound template only: The compound defines the atoms,
       bonds, and geometry, with descriptor tags marking the connection
       atoms::

           path = straight_line(spacing=0.35, N=8, bead_name="PEO")
           monomer = mb.load("C{>}OC{<}", smiles=True)
           compound = path.backmap(templates={"PEO": monomer})

    3. Both: the string owns the bonding rules and the template donates
       only its coordinates (in place of an RDKit embedding); note the
       template needs no tags here. Useful for pairing a relaxed mBuild
       compound for internal geometry with CGsmiles syntax for bonding::

           path = straight_line(spacing=0.35, N=8, bead_name="PEO")
           donor = mb.load("COC", smiles=True)  # e.g. a pre-relaxed conformer
           compound = path.backmap("{#PEO=[>]COC[<]}", templates={"PEO": donor})
    """
    cg_graph = to_cgsmiles_graph(system, fragname_map=fragname_map)
    _, molecule, node_to_beads = _resolve_graph(
        cg_graph,
        fragments,
        legacy=legacy,
        validate=validate,
        templates=templates,
    )

    bead_positions = nx.get_node_attributes(cg_graph, "position")

    # Anchor of each atom: mean position of the bead it descends from
    anchors = {}
    for node in molecule.nodes:
        beads = node_to_beads[node]
        if not beads:
            raise ValueError(
                f"Atom {node} does not trace back to any CG bead; "
                "cannot assign a position."
            )
        anchors[node] = np.mean([bead_positions[b] for b in beads], axis=0)

    if placement not in ("rdkit", "template"):
        raise ValueError(
            f"Argument {placement=} is invalid. "
            "Pass one of 'rdkit', 'template'."
        )
    if placement == "template" and not templates:
        raise ValueError(
            "placement='template' requires the `templates` argument."
        )

    positions = generate_positions(
        molecule=molecule,
        node_to_beads=node_to_beads,
        anchors=anchors,
        bead_fragnames=nx.get_node_attributes(cg_graph, "fragname"),
        placement=placement,
        templates=templates,
        seed=seed,
    )

    return _molecule_to_compound(cg_graph, molecule, node_to_beads, positions)


def _molecule_to_compound(cg_graph, molecule, node_to_beads, positions):
    """Build an mbuild Compound from a resolved molecule graph.

    The hierarchy is one child Compound per CG bead (named after the
    bead) containing that bead's atoms. Atoms descending from several
    beads are placed in the compound of the lowest bead index. If the
    CG bond graph holds multiple disconnected molecules, an
    intermediate ``Molecule`` compound groups each molecule's beads.
    """
    from ele import element_from_symbol
    from ele.exceptions import ElementError

    root = Compound(name="Backmapped")
    bead_compounds = {}
    for node in sorted(cg_graph.nodes):
        bead_compounds[node] = Compound(name=cg_graph.nodes[node]["beadname"])

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

    components = _sorted_components(cg_graph)
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


def _validate_connectivity(cg_graph, molecule, node_to_beads):
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

    expected = {(min(u, v), max(u, v)) for u, v in cg_graph.edges}
    missing = sorted(expected - realized)
    if missing:
        names = nx.get_node_attributes(cg_graph, "beadname")
        details = ", ".join(
            f"({u} '{names[u]}' - {v} '{names[v]}')" for u, v in missing[:10]
        )
        if len(missing) > 10:
            details += ", ..."
        raise ValueError(
            f"{len(missing)} CG bond(s) were not realized as "
            f"atomistic bonds: {details}. This usually means a fragment ran "
            "out of bonding descriptors; each fragment needs at least as many "
            "descriptors as the highest degree of the beads that use it. "
            "Pass validate=False to skip this check."
        )
