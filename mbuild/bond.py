from copy import copy
from warnings import warn
from mbuild.atom import Atom
__author__ = 'sallai'
import numpy as np

class Bond(object):
    @classmethod
    def create(cls, atom1, atom2, kind='undefined', color='black'):
        assert(not atom1 == atom2)
        b = Bond()
        b.kind = kind
        b.color = color
        b.atom1 = atom1
        b.atom2 = atom2
        return b

    def clone(self):
        return copy(self)

    def cloneImage(self):
        b = self.clone()
        b.atom1 = self.atom2
        b.atom2 = self.atom1
        return b

    def cloneWithOrder(self, type_A, type_B):
        ab = self.clone()
        if isinstance(self.atom1, type_A) and isinstance(self.atom2, type_B):
            ab.atom1 = self.atom1
            ab.atom2 = self.atom2
            return ab
        elif isinstance(self.atom1, type_B) and isinstance(self.atom2, type_A):
            ab.atom1 = self.atom2
            ab.atom2 = self.atom1
            return ab
        warn("cannot clone bond " + str(bond) + " with order " + str(type_A) + "," + str(type_B))

    @classmethod
    def orderBond(cls, bond, type_A, type_B):
        cls.cloneWithOrder(bond, type_A, type_B)

    def com(self):
        return [sum(y) / len(y) for y in zip(self.atom1.pos, self.atom2.pos)]

    def hasCommonAtomsWith(self, other):
        assert(isinstance(other, Bond))
        if self.atom1 == other.atom1 or self.atom1 == other.atom2 or self.atom2 == other.atom1 or self.atom2 == other.atom2:
            return True
        else:
            return False

    def hasAtomKinds(self, atomType1, atomType2):
        if isinstance(atomType1, type):
            atomType1 = atomType1.kind
        if isinstance(atomType2, type):
            atomType2 = atomType2.kind

        return (self.atom1.kind == atomType1 and self.atom2.kind == atomType2) or (self.atom1.kind == atomType2 and self.atom2.kind == atomType1)

    def plot(self, ax):
        epsilon = 0.3
        pos1 = np.array(self.atom1.pos)
        pos2 = np.array(self.atom2.pos)
        v12 = pos2 - pos1 # vector from atom1 to atom2
        d12 = np.linalg.norm(v12) # atom1-atom2 distance
        p1 = pos1 + v12/d12 * epsilon
        p2 = pos1 + v12/d12 * (d12 - epsilon)
        ax.plot([p1[0], p2[0]],[p1[1], p2[1]],[p1[2], p2[2]], '-', color=self.color)

    def __hash__(self):
        return hash((self.kind, self.atom1, self.atom2))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()




