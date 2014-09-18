from mbuild.mbase import MBase
from mbuild.part_mixin import PartMixin

__author__ = 'sallai'

from copy import deepcopy

import numpy as np


class Bond(MBase, PartMixin):
    """Connection between two Atoms.

    Attributes:
        atom1 (Atom): First Atom in the bond.
        atom2 (Atom): Second Atom in the bond.
        parent (Compound): Compound to which the Bond belongs.
    """
    __slots__ = ['_atom1', '_atom2', 'kind', 'parent', 'referrers']

    def __init__(self, atom1, atom2, kind=None):
        """Initialize a Bond.

        Args:
            atom1 (Atom): First Atom or Port in the bond.
            atom2 (Atom): Second Atom or Port in the bond.

        """
        super(Bond, self).__init__(atom1, atom2, kind=kind)
        assert(not atom1 == atom2)

        # If a Port is used to initialize a Bond, the Atom that the Port is
        # anchored to will be used to create the Bond.
        from mbuild.port import Port
        if isinstance(atom1, Port):
            atom1 = atom1.anchor
        if isinstance(atom2, Port):
            atom2 = atom2.anchor
        self._atom1 = atom1
        self._atom2 = atom2

        if kind is not None:
            self.kind = kind
        else:
            self.kind = '{0}-{1}'.format(atom1.kind, atom2.kind)

        # self.parent = None

        # Ensure Atoms in Bond know about the Bond.
        atom1.bonds.add(self)
        atom2.bonds.add(self)

        # self.referrers = set()

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    # def ancestors(self):
    #     """Generate all ancestors of the Compound recursively.
    #
    #     Yields:
    #         ancestor (Compound): A Compound one or more levels higher in the
    #             hierarchy.
    #
    #     """
    #     yield self.parent
    #     if self.parent is not None:
    #         for a in self.parent.ancestors():
    #             yield a

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
        return isinstance(bond, Bond) and (self.atom1 == bond.atom1 and self.atom2 == bond.atom2
             or self.atom2 == bond.atom1 and self.atom1 == bond.atom1)

    def __repr__(self):
        return "Bond{0}({1}, {2})".format(id(self), self.atom1, self.atom2)

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        
        # remember the topmost component being deepcopied
        if len(memo) == 0:
            print 'bond is root of deepcopy'
            memo[0] = self
        memo[id(self)] = newone
        
        # copy fields that don't need recursion
        newone.kind = self.kind
        newone.referrers = set()

        # do the rest recursively        
        newone._atom1 = deepcopy(self.atom1, memo)
        newone._atom2 = deepcopy(self.atom2, memo)
        newone._atom1.bonds.add(newone)
        newone._atom2.bonds.add(newone)

        # copy the parent of everybody, except the topmost compound being deepcopied
        if memo[0] == self:
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)
        
        return newone
