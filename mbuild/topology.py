__author__ = 'sallai'
import itertools
from copy import deepcopy

import numpy as np

from mdtraj.core.topology import Topology as MDTTopology
from mdtraj.core.element import Element
from mdtraj.core import element as elem


class Topology(MDTTopology):
    """Derivative of MDTraj's Topology class with additional functionalities.

    Most notably, provides conversion to and from mBuild Compounds.
    """
    def __init__(self):
        """Initialize an mBuild Topology. """
        super(Topology, self).__init__()
        self._ff_bonds = []
        self._ff_angles = []
        self._ff_dihedrals = []
        self._ff_impropers = []

    @property
    def ff_bonds(self):
        return iter(self._ff_bonds)

    @property
    def ff_angles(self):
        return iter(self._ff_angles)

    @property
    def ff_dihedrals(self):
        return iter(self._ff_dihedrals)

    @property
    def ff_impropers(self):
        return iter(self._ff_impropers)

    @property
    def n_ff_bonds(self):
        return sum(1 for _ in self.ff_bonds)

    @property
    def n_ff_angles(self):
        return sum(1 for _ in self.ff_angles)

    @property
    def n_ff_dihedrals(self):
        return sum(1 for _ in self.ff_dihedrals)

    @property
    def n_ff_impropers(self):
        return sum(1 for _ in self.ff_impropers)

    def sort_atoms_alphabetically(self, atoms):
        """Sort a list of atoms alphabetically by their lowercase names. """
        atoms.sort(key=lambda x: x.name.lower())
        return atoms

    def add_ff_bond(self, atom1, atom2):
        atoms = self.sort_atoms_alphabetically([atom1, atom2])
        self._ff_bonds.append(ForcefieldBond(*atoms))

    def add_ff_angle(self, atom1, atom2, atom3):
        atoms = self.sort_atoms_alphabetically([atom1, atom3])
        self._ff_angles.append(ForcefieldAngle(atoms[0], atom2, atoms[1]))

    def add_ff_dihedral(self, atom1, atom2, atom3, atom4):
        atoms = self.sort_atoms_alphabetically([atom1, atom4])
        if atoms[0] == atom1:
            self._ff_dihedrals.append(ForcefieldDihedral(
                atom1, atom2, atom3, atom4))
        elif atoms[1] == atom1:
            self._ff_dihedrals.append(ForcefieldDihedral(
                atom4, atom3, atom2, atom1))

    def load_ff_bonds(self):
        """Convert from (Atom1, Atom2) to ForcefieldBonds. """
        for bond in self.bonds:
            self.add_ff_bond(bond[0], bond[1])

    def enumerate_angles(self, node, neighbors):
        """Find all angles around a node. """
        for pair in itertools.combinations(neighbors, 2):
            self.add_ff_angle(pair[0], node, pair[1])

    def enumerate_dihedrals(self, node_1, neighbors_1, node_2, neighbors_2):
        """Find all dihedrals around a pair of nodes. """
        # We need to make sure we don't remove the node from the neighbor lists
        # that we will be re-using in the following iterations.
        neighbors_1 = set(neighbors_1) - set([node_2])
        neighbors_2.remove(node_1)

        for pair in itertools.product(neighbors_1, neighbors_2):
            self.add_ff_dihedral(pair[0], node_1, node_2, pair[1])

    def find_forcefield_terms(self, bonds=True, angles=True, dihedrals=True,
                              impropers=True):
        """Convert Bonds to ForcefieldBonds and find angles and dihedrals. """
        if bonds:
            self.load_ff_bonds()

        if any([angles, dihedrals, impropers]):
            graph = self.to_bondgraph()
            for node_1 in graph.nodes_iter():
                neighbors_1 = graph.neighbors(node_1)
                if len(neighbors_1) > 1:
                    if angles:
                        self.enumerate_angles(node_1, neighbors_1)
                    if dihedrals:
                        for node_2 in neighbors_1:
                            if node_2.index > node_1.index:
                                neighbors_2 = graph.neighbors(node_2)
                                if len(neighbors_2) > 1:
                                    self.enumerate_dihedrals(
                                        node_1, neighbors_1, node_2, neighbors_2)
                    if impropers:
                        # TODO: implement
                        pass

    @classmethod
    def from_compound(cls, compound, atom_list=None, bond_list=None):
        """Create a Topology from a Compound.

        Args:
            compound:
            atom_list:
            bond_list:
        Returns:
            out (mbuild.Topology):
        """
        out = cls()
        atom_mapping = {}

        c = out.add_chain()
        r = out.add_residue("RES", c)

        if atom_list is None:
            atom_list = compound.atom_list_by_kind('*', excludeG=True)

        for atom in atom_list:
            try:
                e = elem.get_by_symbol(atom.kind)
            except:
                e = Element(1000, atom.kind, atom.kind, 1.0)

            a = out.add_atom(atom.kind, e, r)
            atom_mapping[atom] = a

        if bond_list is None:
            bond_list = compound.bond_list_by_kind(kind='*')

        for idx, bond in enumerate(bond_list):
            a1 = bond.atom1
            a2 = bond.atom2
            out.add_bond(atom_mapping[a1], atom_mapping[a2])
        return out


class ForcefieldBond(object):
    """Connection between two Atoms.

    Attributes:
        atom1 (Atom): First Atom in the bond.
        atom2 (Atom): Second Atom in the bond.
        kind (str, optional): Descriptive name for the bond.
    """

    def __init__(self, atom1, atom2, kind=None):
        """
        """
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

    def other_atom(self, atom):
        """Returns the other Atom in the Bond. """
        if self._atom1 is atom:
            return self._atom2
        elif self._atom2 is atom:
            return self._atom1

    def distance(self, periodicity=np.array([0.0, 0.0, 0.0])):
        """Vectorized distance calculation considering minimum image. """
        d = np.abs(self.atom1 - self.atom2)
        d = np.where(d > 0.5 * periodicity, periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

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

    Attributes:
        atom1 (Atom): First Atom in the angle.
        atom2 (Atom): Second Atom in the angle.
        atom3 (Atom): Third Atom in the angle.
        kind (str, optional): Descriptive name for the angle.
    """

    def __init__(self, atom1, atom2, atom3, kind=None):
        """
        """
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

    Attributes:
        atom1 (Atom): First Atom in the dihedral.
        atom2 (Atom): Second Atom in the dihedral.
        atom3 (Atom): Third Atom in the dihedral.
        atom4 (Atom): Fourth Atom in the dihedral.
        kind (str, optional): Descriptive name for the dihedral.
    """

    def __init__(self, atom1, atom2, atom3, atom4, kind=None):
        """
        """
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
