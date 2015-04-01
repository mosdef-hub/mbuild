import itertools

import simtk.unit as units

from mbuild.orderedset import OrderedSet
from mbuild.tools.parameterize.atomtyper import find_atomtypes


def apply_forcefield(topology, forcefield, debug=True):
    """Apply a forcefield to a Topology. """
    if forcefield.lower() == 'opls-aa':
        from mbuild.tools.parameterize.oplsaa.forcefield import OPLSForcefield
        ff = OPLSForcefield(topology)
    else:
        raise ValueError("Unsupported forcefield: '{0}'".format(forcefield))

    bond_graph = prepare_atoms(topology)
    find_atomtypes(topology._atoms, forcefield, debug=debug)
    ff.resolve_bondtypes()
    enumerate_forcefield_terms(topology, bond_graph)
    topology.gen_pairs(n_excl=4)
    ff = ff.prune()


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
                               dihedrals=True, impropers=True):
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
    atoms.sort(key=lambda x: x.bondtype.lower())
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
            # Sort by inside pair
            if pair[0] == pair[1]:
                atoms = sort_atoms_alphabetically(node_1, node_2)
                if node_1 == atoms[0]:
                    topology.add_ff_dihedral(ForcefieldDihedral(
                        pair[0], node_1, node_2, pair[1]))
                else:
                    topology.add_ff_dihedral(ForcefieldDihedral(
                        pair[1], node_2, node_1, pair[0]))
            # Sort by outside pair
            else:
                atoms = sort_atoms_alphabetically([pair[0], pair[1]])
                if atoms[0] == pair[0]:
                    topology.add_ff_dihedral(ForcefieldDihedral(
                        pair[0], node_1, node_2, pair[1]))
                else:
                    topology.add_ff_dihedral(ForcefieldDihedral(
                        pair[1], node_2, node_1, pair[0]))


def enumerate_impropers(topology, node, neighbors):
    """Find all impropers around a node. """
    for triplet in itertools.combinations(neighbors, 3):
        atoms = sort_atoms_alphabetically([node, triplet[2]])
        if atoms[0] == node:
            topology.add_ff_improper(ForcefieldImproper(
                node, triplet[0], triplet[1], triplet[2]))
        elif atoms[1] == node:
            topology.add_ff_improper(ForcefieldImproper(
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

    def resolve_bondtypes(self):
        """ """
        for atom in self.top.atoms:
            if atom.atomtype not in self.atomtypes:
                print("Could not find atomtype: '{0}' in forcefield.".format(atom.atomtype))
            else:
                atom.bondtype = self.atomtypes[atom.atomtype].bondtype
                if not hasattr(atom, 'charge') or not atom.charge:
                    atom.charge = self.atomtypes[atom.atomtype].charge

    def find_atom_types(self, bondtype):
        """If the id is the atom type, return the AtomType object. """
        matching_atom_types = []
        bondtype = str(bondtype)
        for kind, atomtype in self.atomtypes.items():
            if bondtype.endswith('*'):
                # id is a wildcard ending in *
                prefix = bondtype[:-1]

                if atomtype.bondtype.startswith(prefix):
                    matching_atom_types.append(kind)
                elif kind.startswith(prefix):
                    matching_atom_types.append(kind)
            else:
                # id is not a wildcard
                if bondtype == atomtype.bondtype:
                    matching_atom_types.append(kind)
        return matching_atom_types

    def prune(self):
        """Create force field with only information relevant to topology. """

        bonds_to_remove = set()
        for bond in self.top.ff_bonds:
            if bond.kind not in self.bondtypes:
                print("Could not find bondtype: '{0}' in forcefield.".format(bond.kind))
                bonds_to_remove.add(bond)
        self.top._ff_bonds = self.top._ff_bonds - bonds_to_remove

        angles_to_remove = set()
        for angle in self.top.ff_angles:
            if angle.kind not in self.angletypes:
                print("Could not find angletype: '{0}' in forcefield.".format(angle.kind))
                angles_to_remove.add(angle)
        self.top._ff_angles = self.top._ff_angles - angles_to_remove

        dihedrals_to_remove = set()
        for dihedral in self.top.ff_dihedrals:
            if dihedral.kind not in self.dihedraltypes:
                print("Could not find dihedraltype: '{0}' in forcefield.".format(dihedral.kind))
                dihedrals_to_remove.add(dihedral)
        self.top._ff_dihedrals = self.top._ff_dihedrals - dihedrals_to_remove

        impropers_to_remove = set()
        for improper in self.top.ff_impropers:
            if improper.kind == ('O_2', 'C_2', 'OS', 'CT'):
                print("Keeping explicitcly defined improper: '{0}'".format(improper.kind))
            elif improper.kind not in self.impropertypes:
                print("Could not find impropertype: '{0}' in forcefield.".format(improper.kind))
                impropers_to_remove.add(improper)
        self.top._ff_impropers = self.top._ff_impropers - impropers_to_remove


        # ff = Forcefield(self.topology)
        # ff.atomtypes.update(self.atomtypes)
        #
        # retained_types = [atom.atomtype for atom in self.topology.atoms]
        #
        # # Prune the atom types
        # for atomtype in list(ff.atomtypes.keys()):
        #     if atomtype not in retained_types:
        #         del ff.atomtypes[atomtype]
        #
        #
        # # Prune the bond types, resolving wildcards.
        # for (bondtype1, bondtype2), bond_type in self.bondtypes.items():
        #     atom_types1 = ff.find_atom_types(bondtype1)
        #     atom_types2 = ff.find_atom_types(bondtype2)
        #
        #     # For every combination of matching atom kinds, create a bond type.
        #     for (atom_type1, atom_type2) in itertools.product(atom_types1, atom_types2):
        #         pair = (atom_type1, atom_type2)
        #         ff.bondtypes[pair] = bond_type
        #
        # # Prune the angle types, resolving wildcards.
        # for (bondtype1, bondtype2, bondtype3), angle_type in self.angletypes.items():
        #     atom_types1 = ff.find_atom_types(bondtype1)
        #     atom_types2 = ff.find_atom_types(bondtype2)
        #     atom_types3 = ff.find_atom_types(bondtype3)
        #
        #     # For every combination of matching atom kinds, create an angle type.
        #     for (atom_type1, atom_type2, atom_type3) in itertools.product(atom_types1, atom_types2, atom_types3):
        #         triplet = (atom_type1, atom_type2, atom_type3)
        #         ff.angletypes[triplet] = angle_type
        # return ff


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
    kind : str, optional, default='atom1.bondtype-atom2.bondtype'
        Descriptive name for the bond.

    """
    def __init__(self, atom1, atom2, kind=None):
        self.atom1 = atom1
        self.atom2 = atom2
        if kind:
            self.kind = kind
        else:
            self.kind = (atom1.bondtype, atom2.bondtype)

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
    kind : str, optional, default='atom1.bondtype-atom2.bondtype-atom3.bondtype'
        Descriptive name for the angle.
    """

    def __init__(self, atom1, atom2, atom3, kind=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        if kind:
            self.kind = kind
        else:
            self.kind = (atom1.bondtype, atom2.bondtype, atom3.bondtype)

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
    kind : str, optional, default='atom1.bondtype-atom2.bondtype-atom3.bondtype-atom4.bondtype'
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
            self.kind = (atom1.bondtype, atom2.bondtype, atom3.bondtype, atom4.bondtype)

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
    kind : str, optional, default='atom1.bondtype-atom2.bondtype-atom3.bondtype-atom4.bondtype'
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
            self.kind = (atom1.bondtype, atom2.bondtype, atom3.bondtype, atom4.bondtype)

    def update_kind(self):
        self.kind = (self.atom1.bondtype, self.atom2.bondtype, self.atom3.bondtype, self.atom4.bondtype)

    def __repr__(self):
        return "Dihedral{0}({1}, {2}, {3}, {4})".format(
            id(self), self.atom1, self.atom2, self.atom3, self.atom4)
