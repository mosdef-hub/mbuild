import itertools

import simtk.unit as units

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


class Forcefield(object):
    """ """
    def __init__(self, topology):
        self.topology = self.top = topology
        self.atomtypes = dict()
        self.bondtypes = dict()
        self.angletypes = dict()
        self.dihedraltypes = dict()
        self.impropertypes = dict()

    def find_atom_types(self, atom_id):
        """If the id is the atom type, return the AtomType object. """
        matching_atom_types = []
        atom_id = str(atom_id)
        for kind, atom_type in self.atomtypes.items():
            if atom_id.endswith('*'):
                # id is a wildcard ending in *
                prefix = atom_id[:-1]

                if atom_type.alias.startswith(prefix):
                    matching_atom_types.append(kind)
                elif kind.startswith(prefix):
                    matching_atom_types.append(kind)
            else:
                # id is not a wildcard
                if atom_id == atom_type.alias:
                    matching_atom_types.append(kind)

        return matching_atom_types

    def prune(self):
        """Create force field with only information relevant to topology. """
        ff = Forcefield()
        ff.atomtypes.update(self.atomtypes)

        all_kinds = ff.atomtypes.keys()
        retained_types = [atom.atomtype for atom in self.topology.atoms]

        # Prune the atom types
        for atomtype in all_kinds:
            if atomtype not in retained_types:
                del ff.atomtypes[atomtype]

        # Prune the bond types, resolving wildcards.
        for (alias1, alias2), bond_type in self.bondtypes.items():
            atom_types1 = ff.find_atom_types(alias1)
            atom_types2 = ff.find_atom_types(alias2)

            # For every combination of matching atom kinds, create a bond type.
            for (atom_type1, atom_type2) in itertools.product(atom_types1, atom_types2):
                pair = (atom_type1, atom_type2)
                ff.bondtypes[pair] = bond_type

        # Prune the angle types, resolving wildcards.
        for (alias1, alias2, alias3), angle_type in self.angletypes.items():
            atom_types1 = ff.find_atom_types(alias1)
            atom_types2 = ff.find_atom_types(alias2)
            atom_types3 = ff.find_atom_types(alias3)

            # For every combination of matching atom kinds, create an angle type.
            for (atom_type1, atom_type2, atom_type3) in itertools.product(atom_types1, atom_types2, atom_types3):
                triplet = (atom_type1, atom_type2, atom_type3)
                ff.angletypes[triplet] = angle_type
        return ff


class ForcefieldAtomtype(object):
    """
    """
    def __init__(self, atomtype, bondtype, atomic_number=0,
                 mass=0.0 * units.amu,
                 charge=0.0 * units.elementary_charge,
                 sigma=0.0 * units.angstroms,
                 epsilon=0.0 * units.kilojoules_per_mole):
        self.atomtype = atomtype
        self.bondtype = bondtype
        self.atomic_number = atomic_number
        self.mass = mass
        self.charge = charge
        self.sigma = sigma
        self.epsilon = epsilon


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
        self.atom1 = atom1
        self.atom2 = atom2
        if kind:
            self.kind = kind
        else:
            self.kind = '{0}-{1}'.format(atom1.name, atom2.name)

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
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        if kind:
            self.kind = kind
        else:
            self.kind = '{0}-{1}-{2}'.format(atom1.name, atom2.name, atom3.name)

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
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        if kind:
            self.kind = kind
        else:
            self.kind = '{0}-{1}-{2}-{3}'.format(
                atom1.name, atom2.name, atom3.name, atom4.name)

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
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        if kind:
            self.kind = kind
        else:
            self.kind = '{0}-{1}-{2}-{3}'.format(
                atom1.name, atom2.name, atom3.name, atom4.name)

    def __repr__(self):
        return "Dihedral{0}({1}, {2}, {3}, {4})".format(
            id(self), self.atom1, self.atom2, self.atom3, self.atom4)
