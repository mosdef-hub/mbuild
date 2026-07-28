"""Define CGsmiles fragments from tagged mBuild Compounds."""

import networkx as nx

__all__ = ["compound_to_fragment_graph"]

_DESCRIPTOR_PREFIXES = ("<", ">", "$", "!")


def compound_to_fragment_graph(compound, fragname):
    """Convert a tagged mBuild Compound into a CGsmiles fragment graph.

    This lets an mBuild Compound define a fragment's atomistic
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
        raise ValueError(f"Compound for fragment '{fragname}' has no heavy atoms.")
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
        # mbuild's add_bond stores 0.0 for an unspecified bond order
        if order is None or order == "default" or order == 0:
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
