import itertools

from mbuild.orderedset import OrderedSet

from mbuild.tools.parameterize.atomtyper import find_atomtypes


def apply_forcefield(topology, forcefield, debug=True):
    """Apply a forcefield to a Topology. """

    bond_graph = prepare_atoms(topology)
    find_atomtypes(topology._atoms, forcefield, debug=True)
    enumerate_forcefield_terms(topology, bond_graph)


def prepare_atoms(topology):
    """Add neighbors and white- and blacklists to each atom.

    Note
    ----
    The use of ordered sets is not strictly necessary but it helps when
    debugging because it shows the order in which rules are added.

    Parameters
    ----------
    atoms : list of Atom objects
        The atoms whose atomtypes you are looking for. Atoms must have a
        property `neighbors` which is a list of other atoms that they are
        bonded to.

    """
    bond_graph = topology.to_bondgraph()
    for atom in bond_graph.nodes_iter():
        atom.neighbors = bond_graph.neighbors(atom)
        atom.whitelist = OrderedSet()
        atom.blacklist = OrderedSet()
    return bond_graph


def enumerate_forcefield_terms(topology, graph, bonds=True, angles=True,
                               dihedrals=True, impropers=False):
    """Convert Bonds to ForcefieldBonds and find angles and dihedrals. """
    if bonds:
        convert_to_ff_bonds(topology)

    if any([angles, dihedrals, impropers]):
        for node_1 in graph.nodes_iter():
            neighbors_1 = graph.neighbors(node_1)
            if len(neighbors_1) > 1:
                if angles:
                    enumerate_angles(topology, node_1, neighbors_1)
                if dihedrals:
                    for node_2 in neighbors_1:
                        if node_2.index > node_1.index:
                            neighbors_2 = graph.neighbors(node_2)
                            if len(neighbors_2) > 1:
                                enumerate_dihedrals(topology,
                                    node_1, neighbors_1, node_2, neighbors_2)
                if impropers and len(neighbors_1) >= 3:
                    enumerate_impropers(topology, node_1, neighbors_1)


def sort_atoms_alphabetically(atoms):
    """Sort a list of atoms alphabetically by their lowercase names. """
    atoms.sort(key=lambda x: x.name.lower())
    return atoms


def convert_to_ff_bonds(topology):
    """Convert from tuples of (Atom1, Atom2) to ForcefieldBonds. """
    for bond in topology.bonds:
        atoms = sort_atoms_alphabetically([bond[0], bond[1]])
        topology.add_ff_bond(ForcefieldBond(*atoms))


def enumerate_angles(topology, node, neighbors):
    """Find all angles around a node. """
    for pair in itertools.combinations(neighbors, 2):
        atoms = sort_atoms_alphabetically([pair[0], pair[1]])
        topology.add_ff_angle(ForcefieldAngle(atoms[0], node, atoms[1]))


def enumerate_dihedrals(topology, node_1, neighbors_1, node_2, neighbors_2):
    """Find all dihedrals around a pair of nodes. """
    # We need to make sure we don't remove the node from the neighbor lists
    # that we will be re-using in the following iterations.
    neighbors_1 = set(neighbors_1) - {node_2}
    neighbors_2.remove(node_1)

    for pair in itertools.product(neighbors_1, neighbors_2):
        if pair[0] != pair[1]:
            atoms = sort_atoms_alphabetically([pair[0], pair[1]])
            if atoms[0] == pair[0]:
                topology.add_ff_dihedral(ForcefieldDihedral(
                    pair[0], node_1, node_2, pair[1]))
            elif atoms[1] == pair[1]:
                topology.add_ff_dihedral.append(ForcefieldDihedral(
                    pair[1], node_2, node_1, pair[0]))


def enumerate_impropers(topology, node, neighbors):
    """Find all impropers around a node. """
    for triplet in itertools.combinations(neighbors, 3):
        atoms = sort_atoms_alphabetically([node, triplet[2]])
        if atoms[0] == node:
            topology.add_ff_improper(ForcefieldImproper(
                node, triplet[0], triplet[1], triplet[2]))
        elif atoms[1] == node:
            topology._ff_dihedrals.append(ForcefieldImproper(
                triplet[2], triplet[1], triplet[0], node))



class ForcefieldBond(object):
    """Connection between two Atoms.

    Attributes
    ----------
    atom1 : mb.Atom
        First Atom in the bond.
    atom2 : mb.Atom
        Second Atom in the bond.
    kind : str, optional, default='atom1.name-atom2.name'
        Descriptive name for the bond.

    """
    def __init__(self, atom1, atom2, kind=None):
        self._atom1 = atom1
        self._atom2 = atom2
        if kind:
            self.kind = kind
        else:
            self.kind = '{0}-{1}'.format(atom1.name, atom2.name)

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    def __hash__(self):
        return id(self.atom1) ^ id(self.atom2)

    def __eq__(self, bond):
        return (isinstance(bond, ForcefieldBond) and
                (self.atom1 == bond.atom1 and self.atom2 == bond.atom2 or
                 self.atom2 == bond.atom1 and self.atom1 == bond.atom1))

    def __repr__(self):
        return "Bond{0}({1}, {2})".format(id(self), self.atom1, self.atom2)


class ForcefieldAngle(object):
    """Triplet formed by three Atoms and two consecutive Bonds.

    Attributes
    ----------
    atom1 : mb.Atom
        First Atom in the angle.
    atom2 : mb.Atom
        Second Atom in the angle.
    atom3 : mb.Atom
        Second Atom in the angle.
    kind : str, optional, default='atom1.name-atom2.name-atom3.name'
        Descriptive name for the angle.
    """

    def __init__(self, atom1, atom2, atom3, kind=None):
        self._atom1 = atom1
        self._atom2 = atom2
        self._atom3 = atom3
        if kind:
            self.kind = kind
        else:
            self.kind = '{0}-{1}-{2}'.format(atom1.name, atom2.name, atom3.name)

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def atom3(self):
        return self._atom3

    def __hash__(self):
        return id(self.atom1) ^ id(self.atom2) ^ id(self.atom3)

    def __eq__(self, angle):
        return (isinstance(angle, ForcefieldAngle) and
                (self.atom1 == angle.atom1 and self.atom2 == angle.atom2 and self.atom3 == angle.atom3 or
                 self.atom3 == angle.atom1 and self.atom2 == angle.atom2 and self.atom1 == angle.atom3))

    def __repr__(self):
        return "Angle{0}({1}, {2}, {3})".format(
            id(self), self.atom1, self.atom2, self.atom3)


class ForcefieldDihedral(object):
    """Quadruplet formed by four Atoms, three Bonds and two angles.

    Attributes
    ----------
    atom1 : mb.Atom
        First Atom in the dihedral.
    atom2 : mb.Atom
        Second Atom in the dihedral.
    atom3 : mb.Atom
        Second Atom in the dihedral.
    atom4 : mb.Atom
        Second Atom in the dihedral.
    kind : str, optional, default='atom1.name-atom2.name-atom3.name-atom4.name'
        Descriptive name for the dihedral.

    """

    def __init__(self, atom1, atom2, atom3, atom4, kind=None):
        self._atom1 = atom1
        self._atom2 = atom2
        self._atom3 = atom3
        self._atom4 = atom4
        if kind:
            self.kind = kind
        else:
            self.kind = '{0}-{1}-{2}-{3}'.format(
                atom1.name, atom2.name, atom3.name, atom4.name)

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def atom3(self):
        return self._atom3

    @property
    def atom4(self):
        return self._atom4

    def __hash__(self):
        return id(self.atom1) ^ id(self.atom2) ^ id(self.atom3) ^ id(self.atom4)

    def __eq__(self, dihedral):
        return (isinstance(dihedral, ForcefieldDihedral) and
                (self.atom1 == dihedral.atom1 and self.atom2 == dihedral.atom2 and self.atom3 == dihedral.atom3 and self.atom4 == dihedral.atom4 or
                 self.atom4 == dihedral.atom1 and self.atom3 == dihedral.atom2 and self.atom2 == dihedral.atom3 and self.atom1 == dihedral.atom4))

    def __repr__(self):
        return "Dihedral{0}({1}, {2}, {3}, {4})".format(
            id(self), self.atom1, self.atom2, self.atom3, self.atom4)


class ForcefieldImproper(object):
    """Quadruplet formed by four Atoms, three Bonds and two angles.

    Attributes
    ----------
    atom1 : mb.Atom
        First Atom in the dihedral.
    atom2 : mb.Atom
        Second Atom in the dihedral.
    atom3 : mb.Atom
        Second Atom in the dihedral.
    atom4 : mb.Atom
        Second Atom in the dihedral.
    kind : str, optional, default='atom1.name-atom2.name-atom3.name-atom4.name'
        Descriptive name for the dihedral.

    """

    def __init__(self, atom1, atom2, atom3, atom4, kind=None):
        self._atom1 = atom1
        self._atom2 = atom2
        self._atom3 = atom3
        self._atom4 = atom4
        if kind:
            self.kind = kind
        else:
            self.kind = '{0}-{1}-{2}-{3}'.format(
                atom1.name, atom2.name, atom3.name, atom4.name)

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def atom3(self):
        return self._atom3

    @property
    def atom4(self):
        return self._atom4

    def __hash__(self):
        return id(self.atom1) ^ id(self.atom2) ^ id(self.atom3) ^ id(self.atom4)

    def __eq__(self, dihedral):
        return (isinstance(dihedral, ForcefieldImproper) and
                (self.atom1 == dihedral.atom1 and self.atom2 == dihedral.atom2 and self.atom3 == dihedral.atom3 and self.atom4 == dihedral.atom4 or
                 self.atom4 == dihedral.atom1 and self.atom3 == dihedral.atom2 and self.atom2 == dihedral.atom3 and self.atom1 == dihedral.atom4))

    def __repr__(self):
        return "Dihedral{0}({1}, {2}, {3}, {4})".format(
            id(self), self.atom1, self.atom2, self.atom3, self.atom4)
