"""Forward coarse-graining: atomistic Compounds to CG Paths.

This is the inverse of ``backmap``: an atomistic ``mbuild.Compound`` is
partitioned into CG beads and returned as an ``mbuild.path.Path`` (bead
names, bead coordinates, CG bond graph) together with an explicit
atom-to-bead mapping. The same CGsmiles fragment string that drives
backmapping can drive the forward direction, so one string describes
the resolution transform both ways.

Three ways to define the mapping:

1. CGsmiles fragments (``fragments=...`` and/or ``templates=...``):
   the compound's bond graph is tiled with the fragment graphs by
   subgraph matching. Reproducible and declarative; the headline
   route.
2. Hierarchy (``beads=["PEO", ...]``): sub-compounds with a listed
   name become beads (with all their particles). This is how compounds
   produced by ``backmap``, which have one sub-compound per bead,
   coarse-grain back trivially.
3. Explicit (``mapping={particle: bead_index}``).


This process utilizes the CGSmiles python package.
If you use this feature please cite:

CGsmiles: A Versatile Line Notation for Molecular Representations
across Multiple Resolutions. Fabian Grünewald, Leif Seute, Riccardo
Alessandri, Melanie König, and Peter C. Kroon. Journal of Chemical
Information and Modeling 2025 65 (7), 3405-3419.
DOI: 10.1021/acs.jcim.5c00064
"""

from collections import defaultdict
from dataclasses import dataclass, field

import networkx as nx
import numpy as np

from .convert import _sorted_components
from .fragments import _particle_symbol, compound_to_fragment_graph

__all__ = ["coarse_grain", "CGMapping"]


@dataclass
class CGMapping:
    """Correspondence between an atomistic Compound and its CG beads.

    Atoms are identified by their index in ``list(compound.particles())``
    (the order mBuild iterates particles in); beads by their node index
    in the returned Path. This indexed form is directly applicable to
    trajectory frames written from the same compound.

    Attributes
    ----------
    atom_to_bead : dict[int, int]
        Particle index -> bead index.
    bead_to_atoms : dict[int, list[int]]
        Bead index -> sorted particle indices belonging to it.
    """

    atom_to_bead: dict = field(default_factory=dict)
    bead_to_atoms: dict = field(default_factory=dict)

    @property
    def n_beads(self):
        """Number of CG beads."""
        return len(self.bead_to_atoms)


def coarse_grain(
    compound,
    fragments=None,
    templates=None,
    beads=None,
    mapping=None,
    bead_names=None,
    center="geometry",
):
    """Coarse-grain an atomistic Compound into a CG Path.

    Exactly one mapping mode must be given: ``fragments`` (and/or
    ``templates``), ``beads``, or ``mapping``.

    Parameters
    ----------
    compound : mbuild.Compound, required
        The atomistic compound to coarse-grain.
    fragments : str, optional
        A CGsmiles fragment string, e.g. ``"{#PEO=[>]COC[<]}"``; the
        same string used for backmapping. Each fragment's heavy-atom
        graph is matched into the compound's bond graph (elements,
        charges, and bond orders must agree; hydrogens are folded into
        each heavy atom's hydrogen count) and the compound must tile
        exactly into fragment instances. Hydrogen deficits relative to
        the saturated fragment must be accounted for by inter-fragment
        bonds, and each fragment atom needs enough bonding descriptors
        for the inter-fragment bonds it forms. Bead names are the
        fragment names. Multi-resolution strings are not supported.
    templates : dict[str, mbuild.Compound], optional
        Tagged compounds defining fragments not present in
        ``fragments`` (see ``compound_to_fragment_graph``); usable with
        or without a fragment string, exactly as in ``backmap``.
    beads : list of str, optional
        Hierarchy mode: walking the compound tree from the top, any
        sub-compound whose name is in this list becomes one bead
        containing all its particles (no descent into it). Every
        particle must end up in a bead. Bead names are the sub-compound
        names. Compounds produced by ``backmap`` coarse-grain back with
        their bead names here.
    mapping : dict[mbuild.Compound, int], optional
        Explicit mode: maps every particle of ``compound`` to a bead
        index. Bead indices must be ``0..n_beads-1``.
    bead_names : array-like of str, optional
        Bead names for explicit ``mapping`` mode (one per bead).
        Defaults to ``"_A"``. Ignored in the other modes.
    center : str, default "geometry"
        ``"geometry"`` places each bead at the arithmetic mean of its
        atoms' positions; ``"mass"`` uses the mass-weighted mean (every
        particle then needs a mass).

    Returns
    -------
    path : mbuild.path.Path
        The CG path: one node per bead with bead names and coordinates,
        and a CG bond between beads connected by at least one
        atomistic bond.
    mapping : CGMapping
        Atom/bead correspondence by particle index.

    Notes
    -----
    Bead coordinates are a snapshot computed from the current particle
    positions; re-run ``coarse_grain`` (or apply the returned mapping
    to a trajectory) after the compound moves.

    Examples
    --------
    1. CGsmiles fragment string: the inverse of backmap. Each fragment's
       heavy-atom graph is matched into the compound, which is tiled into
       beads named after the fragments::

           path = straight_line(spacing=0.35, N=8, bead_name="PEO")
           compound = path.backmap("{#PEO=[>]COC[<]}")
           cg_path, mapping = compound.coarse_grain("{#PEO=[>]COC[<]}")

    2. mBuild compound template: a tagged compound defines the fragment
       in place of a fragment string, exactly as in backmap::

           template = mb.load("C{>}OC{<}", smiles=True)
           cg_path, mapping = compound.coarse_grain(templates={"PEO": template})

    3. Hierarchy mode: any sub-compound whose name is listed becomes one
       bead. Compounds produced by backmap carry their bead names, so
       they coarse-grain straight back::

           cg_path, mapping = compound.coarse_grain(beads=["PEO"])
    """
    from mbuild.path.build import Path
    from mbuild.path.namers import BEAD_NAME_DTYPE

    modes = [fragments is not None or templates, beads is not None, mapping is not None]
    if sum(bool(m) for m in modes) != 1:
        raise ValueError(
            "Provide exactly one mapping mode: `fragments`/`templates`, "
            "`beads`, or `mapping`."
        )
    if center not in ("geometry", "mass"):
        raise ValueError(f"Argument {center=} is invalid. Pass 'geometry' or 'mass'.")

    particles = list(compound.particles())
    index_of = {particle: index for index, particle in enumerate(particles)}

    if mapping is not None:
        bead_members, names = _beads_from_explicit(
            particles, index_of, mapping, bead_names
        )
    elif beads is not None:
        bead_members, names = _beads_from_hierarchy(compound, index_of, beads)
    else:
        bead_members, names = _beads_from_fragments(
            compound, particles, index_of, fragments, templates
        )

    atom_to_bead = {}
    for bead_index, members in enumerate(bead_members):
        for atom in members:
            atom_to_bead[atom] = bead_index

    coordinates = np.array(
        [
            _bead_center([particles[i] for i in members], center)
            for members in bead_members
        ]
    )

    cg_graph = nx.Graph()
    cg_graph.add_nodes_from(range(len(bead_members)))
    for particle_a, particle_b in compound.bonds():
        bead_a = atom_to_bead[index_of[particle_a]]
        bead_b = atom_to_bead[index_of[particle_b]]
        if bead_a != bead_b:
            cg_graph.add_edge(bead_a, bead_b)

    path = Path(
        coordinates=coordinates,
        bond_graph=cg_graph,
        bead_name=np.array(names, dtype=BEAD_NAME_DTYPE),
    )
    cg_mapping = CGMapping(
        atom_to_bead=atom_to_bead,
        bead_to_atoms={
            bead: sorted(members) for bead, members in enumerate(bead_members)
        },
    )
    return path, cg_mapping


def _bead_center(members, center):
    """Coordinate of one bead: the (mass-weighted) centroid of its particles.

    With ``center="geometry"`` this is the arithmetic mean of the member
    positions; with ``center="mass"`` the mass-weighted mean, requiring a
    positive total mass over the bead.
    """
    positions = np.array([particle.pos for particle in members], dtype=float)
    if center == "geometry":
        return positions.mean(axis=0)
    masses = []
    for particle in members:
        if particle.mass is None:
            raise ValueError(
                f"center='mass' requires a mass on every particle; "
                f"'{particle.name}' has none. Use center='geometry' instead."
            )
        masses.append(particle.mass)
    masses = np.asarray(masses, dtype=float)
    if masses.sum() <= 0:
        raise ValueError(
            "center='mass' requires a positive total mass per bead; got "
            f"{masses.sum()} for a bead of {len(members)} particle(s). "
            "Use center='geometry' instead."
        )
    return np.average(positions, axis=0, weights=masses)


def _beads_from_explicit(particles, index_of, mapping, bead_names):
    """Group particle indices into beads from an explicit particle to
    bead-index mapping, checking that it covers every particle and uses
    contiguous indices. Returns the per-bead member lists and names.
    """
    missing = [p.name for p in particles if p not in mapping]
    if missing:
        raise ValueError(
            f"`mapping` must cover every particle; {len(missing)} particle(s) "
            f"are missing (first few: {missing[:5]})."
        )
    bead_indices = sorted({int(b) for b in mapping.values()})
    if bead_indices != list(range(len(bead_indices))):
        raise ValueError(
            "`mapping` bead indices must be contiguous starting at 0; "
            f"got {bead_indices[:10]}."
        )
    members = [[] for _ in bead_indices]
    for particle, bead in mapping.items():
        if particle in index_of:
            members[int(bead)].append(index_of[particle])
    if bead_names is None:
        names = ["_A"] * len(members)
    else:
        names = [str(name) for name in bead_names]
        if len(names) != len(members):
            raise ValueError(
                f"`bead_names` has {len(names)} entries for {len(members)} beads."
            )
    return members, names


def _beads_from_hierarchy(compound, index_of, beads):
    """Group particles into beads by walking the compound tree, treating
    each sub-compound whose name is listed in ``beads`` as one bead.
    Returns the per-bead member lists and names, raising if any particle
    is left unassigned.
    """
    bead_set = set(beads)
    bead_compounds = []

    def collect(comp):
        """Recursively descend the tree, collecting sub-compounds whose
        name marks them as a bead.
        """
        for child in comp.children or []:
            if child.name in bead_set:
                bead_compounds.append(child)
            elif child.children:
                collect(child)

    if compound.name in bead_set:
        bead_compounds.append(compound)
    else:
        collect(compound)

    members = []
    names = []
    seen = set()
    for bead_compound in bead_compounds:
        indices = [index_of[p] for p in bead_compound.particles()]
        members.append(indices)
        names.append(bead_compound.name)
        seen.update(indices)
    unmapped = len(index_of) - len(seen)
    if unmapped:
        raise ValueError(
            f"{unmapped} particle(s) are not contained in any sub-compound "
            f"named in `beads` ({sorted(bead_set)}). Every particle must "
            "belong to a bead."
        )
    return members, names


def _beads_from_fragments(compound, particles, index_of, fragments, templates):
    """Tile the compound's heavy-atom graph into fragment instances,
    folding each heavy atom's hydrogens back in, to make one bead per
    instance. Returns the per-bead member lists and fragment names.
    """
    fragment_dict = _fragment_dict(fragments, templates)
    heavy_graph, h_of = _heavy_graph(compound, particles, index_of)

    members = []
    names = []
    for component in _sorted_components(heavy_graph):
        instances = _cover_component(heavy_graph.subgraph(component), fragment_dict)
        for fragname, atoms in instances:
            atom_indices = sorted(atoms)
            for atom in atoms:
                atom_indices.extend(h_of[atom])
            members.append(sorted(atom_indices))
            names.append(fragname)

    # H atoms not bonded to any heavy atom cannot be assigned
    assigned = {i for m in members for i in m}
    unassigned = [particles[i].name for i in range(len(particles)) if i not in assigned]
    if unassigned:
        raise ValueError(
            f"{len(unassigned)} particle(s) could not be assigned to a bead "
            f"(first few: {unassigned[:5]}); hydrogens must be bonded to a "
            "heavy atom."
        )
    return members, names


def _fragment_dict(fragments, templates):
    """Build the fragment-name to fragment-graph dict used for matching,
    from a single all-atom CGsmiles string and/or tagged template
    compounds. Multi-resolution strings are rejected.
    """
    if fragments is not None:
        import re

        from cgsmiles.resolve import MoleculeResolver

        if not isinstance(fragments, str):
            raise NotImplementedError(
                "Forward coarse-graining supports a single all-atom "
                "fragment string; multi-resolution strings are not "
                "supported."
            )
        elements = re.findall(r"\{[^\}]+\}", fragments)
        fragment_dicts = MoleculeResolver.read_fragment_strings(
            elements, last_all_atom=True
        )
        if len(fragment_dicts) != 1:
            raise NotImplementedError(
                "Forward coarse-graining supports a single all-atom "
                "fragment string; multi-resolution strings are not "
                "supported."
            )
        fragment_dict = dict(fragment_dicts[-1])
    else:
        fragment_dict = {}
    if templates:
        for fragname, template in templates.items():
            if fragname not in fragment_dict:
                fragment_dict[fragname] = compound_to_fragment_graph(template, fragname)
    return fragment_dict


def _heavy_graph(compound, particles, index_of):
    """Heavy-atom graph of the compound with hydrogens folded in.

    Nodes are particle indices with ``element``, ``charge``, and
    ``hcount`` (actual bonded hydrogens); edges carry the bond
    ``order``. Also returns a dict mapping each heavy atom to the
    indices of its bonded hydrogens.
    """
    graph = nx.Graph()
    heavy = [
        index
        for index, particle in enumerate(particles)
        if _particle_symbol(particle) != "H"
    ]
    if not heavy:
        raise ValueError("Compound has no heavy atoms to coarse-grain.")
    for index in heavy:
        particle = particles[index]
        graph.add_node(
            index,
            element=_particle_symbol(particle),
            charge=int(particle.charge or 0),
            hcount=0,
        )
    h_of = defaultdict(list)
    for particle_a, particle_b, data in compound.bonds(return_bond_order=True):
        a, b = index_of[particle_a], index_of[particle_b]
        a_is_h = _particle_symbol(particle_a) == "H"
        b_is_h = _particle_symbol(particle_b) == "H"
        if a_is_h or b_is_h:
            if not a_is_h:
                graph.nodes[a]["hcount"] += 1
                h_of[a].append(b)
            elif not b_is_h:
                graph.nodes[b]["hcount"] += 1
                h_of[b].append(a)
            continue
        order = data.get("bond_order", data.get("bo"))
        # mbuild's add_bond stores 0.0 for an unspecified bond order
        if order is None or order == "default" or order == 0:
            order = 1
        graph.add_edge(a, b, order=float(order))
    return graph, h_of


def _frag_order(data):
    """Return a fragment edge's bond order as a float, defaulting to a
    single bond.
    """
    order = data.get("order", 1)
    return float(order)


def _cover_component(component_graph, fragment_dict):
    """Tile one connected heavy-atom component with fragment instances.

    Enumerates node-induced subgraph matches of every fragment, keeping
    only locally consistent ones: in any exact cover, an atom's
    inter-fragment bond count is exactly its heavy degree minus its
    degree inside the instance, so the hydrogen deficit relative to the
    saturated fragment must equal that number and stay within the
    atom's bonding-descriptor capacity. A backtracking search then
    builds an exact cover of the component's atoms, pruning placements
    whose bonds to already-placed instances lack a compatible bonding
    descriptor pair (``<x`` bonds ``>x``, ``$x`` bonds ``$x``); this
    is what rejects frame-shifted tilings that are otherwise
    self-consistent. Descriptor pairs are checked per bond, not
    assigned injectively across an atom's bonds. Returns a list of
    ``(fragname, atom_set)`` sorted by lowest atom index.
    """
    from networkx.algorithms.isomorphism import GraphMatcher

    def node_match(mol_node, frag_node):
        """Atom comparison for the subgraph matcher: element and charge
        must match and the molecule atom's hydrogen count must not exceed
        the fragment atom's.
        """
        return (
            mol_node["element"] == frag_node.get("element")
            and mol_node["charge"] == frag_node.get("charge", 0)
            and mol_node["hcount"] <= frag_node.get("hcount", 0)
        )

    def edge_match(mol_edge, frag_edge):
        """Bond comparison for the subgraph matcher: bond orders must be
        equal.
        """
        return float(mol_edge.get("order", 1)) == _frag_order(frag_edge)

    degree = dict(component_graph.degree)

    def locally_consistent(fragment, match):
        """Keep only matches where every atom's hydrogen deficit equals
        its number of inter-fragment bonds and stays within its
        bonding-descriptor capacity.
        """
        internal = component_graph.subgraph(match)
        for atom, frag_node in match.items():
            external = degree[atom] - internal.degree(atom)
            deficit = (
                fragment.nodes[frag_node].get("hcount", 0)
                - component_graph.nodes[atom]["hcount"]
            )
            if deficit != external:
                return False
            if external > len(fragment.nodes[frag_node].get("bonding", [])):
                return False
        return True

    # Enumerate candidate placements (deduplicated: automorphic matches
    # with identical per-atom descriptor lists are equivalent)
    embeddings = []
    seen = set()
    for fragname, fragment in fragment_dict.items():
        matcher = GraphMatcher(
            component_graph, fragment, node_match=node_match, edge_match=edge_match
        )
        for match in matcher.subgraph_isomorphisms_iter():
            # match: component atom -> fragment node
            if not locally_consistent(fragment, match):
                continue
            profile = tuple(
                sorted(
                    (atom, tuple(fragment.nodes[frag_node].get("bonding", [])))
                    for atom, frag_node in match.items()
                )
            )
            key = (fragname, profile)
            if key in seen:
                continue
            seen.add(key)
            descriptors = {
                atom: fragment.nodes[frag_node].get("bonding", [])
                for atom, frag_node in match.items()
            }
            embeddings.append((fragname, frozenset(match), descriptors))

    atoms = frozenset(component_graph.nodes)
    by_atom = defaultdict(list)
    for embedding_index, (_, emb_atoms, _) in enumerate(embeddings):
        for atom in emb_atoms:
            by_atom[atom].append(embedding_index)

    uncoverable = sorted(atom for atom in atoms if not by_atom[atom])
    if uncoverable:
        details = ", ".join(
            f"{atom} ({component_graph.nodes[atom]['element']}, "
            f"{component_graph.nodes[atom]['hcount']} H)"
            for atom in uncoverable[:5]
        )
        raise ValueError(
            f"No fragment in {sorted(fragment_dict)} consistently matches "
            f"heavy atom(s) {details}. Check elements, bond orders, and "
            "hydrogen counts against the fragment SMILES, and that each "
            "fragment atom has a bonding descriptor per inter-fragment "
            "bond it must form."
        )

    def compatible_with_placed(emb_atoms, descriptors, atom_descriptors):
        """Descriptor-pair check on bonds into already-placed instances."""
        for atom in emb_atoms:
            for neighbor in component_graph.neighbors(atom):
                if neighbor in emb_atoms or neighbor not in atom_descriptors:
                    continue
                if not any(
                    _descriptors_compatible(descriptor_a, descriptor_b)
                    for descriptor_a in descriptors[atom]
                    for descriptor_b in atom_descriptors[neighbor]
                ):
                    return False
        return True

    def search(uncovered, atom_descriptors):
        """Backtracking search that builds an exact cover of the
        component's atoms from the candidate placements, respecting
        bonding-descriptor compatibility. Returns the chosen embedding
        indices, or None if no cover exists.
        """
        if not uncovered:
            return []
        target = min(uncovered)
        for embedding_index in by_atom[target]:
            fragname, emb_atoms, descriptors = embeddings[embedding_index]
            if not uncovered.issuperset(emb_atoms):
                continue
            if not compatible_with_placed(emb_atoms, descriptors, atom_descriptors):
                continue
            atom_descriptors.update(descriptors)
            rest = search(uncovered - emb_atoms, atom_descriptors)
            if rest is not None:
                return [embedding_index] + rest
            for atom in emb_atoms:
                del atom_descriptors[atom]
        return None

    cover = search(atoms, {})
    if cover is None:
        raise ValueError(
            f"Could not tile a molecule of {len(atoms)} heavy atoms with the "
            f"fragments {sorted(fragment_dict)}: no exact cover satisfies "
            "the hydrogen accounting and bonding descriptor compatibility. "
            "Each inter-fragment bond must replace one hydrogen on each "
            "side of the saturated fragments and join a compatible "
            "descriptor pair."
        )
    return sorted(
        ((embeddings[i][0], set(embeddings[i][1])) for i in cover),
        key=lambda instance: min(instance[1]),
    )


def _descriptors_compatible(descriptor_a, descriptor_b):
    """BigSmiles-style descriptor compatibility.

    Strips the bond-order suffix; ``$x`` bonds ``$x``, ``<x`` bonds
    ``>x`` (labels must match exactly).
    """
    a = descriptor_a.rstrip("0123456789")
    b = descriptor_b.rstrip("0123456789")
    if a.startswith("$"):
        return a == b
    if a.startswith("<"):
        return b == ">" + a[1:]
    if a.startswith(">"):
        return b == "<" + a[1:]
    return False
