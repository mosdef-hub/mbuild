__author__ = 'sallai'
import numpy as np

from mdtraj.core.element import Element
from mdtraj.core.topology import Topology as MDTTopology
from mdtraj.core import element as elem


class Topology(MDTTopology):

    def __init__(self):
        super(Topology, self).__init__()
        self._ff_bonds = []
        self._ff_angles = []
        self._ff_dihedrals = []


    @property
    def ff_bonds(self):
        return iter(self._ff_bonds)

    @property
    def ff_angles(self):
        return iter(self._ff_angles)

    @property
    def ff_dihedrals(self):
        return iter(self._ff_dihedrals)

    def add_ff_bond(self, atom1, atom2):
        print "Adding bond: {}-{}".format(atom1, atom2)
        self._ff_bonds.append(ForcefieldBond(atom1, atom2))

    def add_angle(self, atom1, atom2, atom3):
        print "Adding angle: {}-{}-{}".format(atom1, atom2, atom3)
        self._ff_angles.append(ForcefieldAngle(atom1, atom2, atom3))

    def add_dihedral(self, atom1, atom2, atom3, atom4):
        print "Adding dihedral: {}-{}-{}".format(atom1, atom2, atom3, atom4)
        self._ff_dihedrals.append(ForcefieldDihedral(atom1, atom2, atom3, atom4))

    @classmethod
    def from_compound(cls, compound, atom_list=None, bond_list=None):

        out = cls()
        atom_mapping = {}

        c = out.add_chain()
        r = out.add_residue("RES", c)

        if atom_list is None:
            atom_list, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

        for atom in atom_list:
            e = None
            try:
                e = elem.get_by_symbol(atom.kind)
            except:
                e = Element(  1000, atom.kind, atom.kind, 1.0)

            a = out.add_atom(atom.kind, e, r)
            atom_mapping[atom] = a

        if bond_list is None:
            bond_list, bond_id_to_idx = compound.bond_list_by_kind('*', with_id_to_idx_mapping=True)
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