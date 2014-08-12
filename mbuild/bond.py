
__author__ = 'sallai'

from copy import copy, deepcopy
from warnings import warn
import numpy as np
class Bond(object):

    __slots__ = ['_atom1', '_atom2', 'parent', 'referrers']

    def __init__(self, atom1, atom2):
        assert(not atom1 == atom2)
        from mbuild.port import Port

        if isinstance(atom1, Port):
            atom1 = atom1.anchor
        if isinstance(atom2, Port):
            atom2 = atom2.anchor

        self._atom1 = atom1
        self._atom2 = atom2

        self.parent = None

        atom1.bonds.add(self)
        atom2.bonds.add(self)

    @property
    def atom1(self):
        return self._atom1

    @atom1.setter
    def atom1(self, v):
        raise TypeError

    @property
    def atom2(self):
        return self._atom2

    @atom2.setter
    def atom2(self, v):
        raise TypeError

    def ancestors(self):
        """
        Generate all ancestors of the Compound recursively
        :yield ancestors
        """
        yield self.parent
        if self.parent is not None:
            for a in self.parent.ancestors():
                yield a

    def __hash__(self):
        return id(self.atom1) ^ id(self.atom2) 

    def distance(self, periodicity=np.array([0.0, 0.0, 0.0])):
        """Vectorized distance calculation considering minimum image
        """
        d = np.abs(self.atom1 - self.atom2)
        d = np.where(d > 0.5 * periodicity, periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

    def __eq__(self, bond):
        return isinstance(bond, Bond) and (self.atom1 == bond.atom1 and self.atom2 == bond.atom2
             or self.atom2 == bond.atom1 and self.atom1 == bond.atom1)

    def __repr__(self):
        return "Bond{0}({1}, {2})".format(id(self), self.atom1, self.atom2)

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        if len(memo) == 0:
            print 'bond is root of deepcopy'
            memo[0] = self
        memo[id(self)] = newone

        newone._atom1 = deepcopy(self.atom1, memo)
        newone._atom2 = deepcopy(self.atom2, memo)
        newone._atom1.bonds.add(newone)
        newone._atom2.bonds.add(newone)

        # copy the parent of everybody, except the topmost compound being deepcopied
        if memo[0] == self:
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)

        newone.referrers = set()

        return newone

