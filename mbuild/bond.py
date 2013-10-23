__author__ = 'sallai'
import numpy as np

class Bond(object):
    @classmethod
    def create(cls, atom1, atom2, bondType='undefined', color='black'):
        b = Bond()
        b.bondType = bondType
        b.color = color
        # atom1 is the one with the lower hash
        if atom1.__hash__() < atom2.__hash__():
            b.atom1 = atom1
            b.atom2 = atom2
        else:
            b.atom1 = atom2
            b.atom2 = atom1

        return b

    def hasCommonAtomsWith(self, other):
        assert(isinstance(other, Bond))
        if self.atom1 == other.atom1 or self.atom1 == other.atom2 or self.atom2 == other.atom1 or self.atom2 == other.atom2:
            return True
        else:
            return False

    def __hash__(self):
        # atom1 and atom2 are interchangeable
        return (self.atom1.__hash__() * self.atom2.__hash__()) ^ self.atom1.__hash__() ^ self.atom2.__hash__()

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def plot(self, ax):
        epsilon = 0.3
        pos1 = np.array(self.atom1.pos)
        pos2 = np.array(self.atom2.pos)
        v12 = pos2 - pos1 # vector from atom1 to atom2
        d12 = np.linalg.norm(v12) # atom1-atom2 distance
        p1 = pos1 + v12/d12 * epsilon
        p2 = pos1 + v12/d12 * (d12 - epsilon)
        ax.plot([p1[0], p2[0]],[p1[1], p2[1]],[p1[2], p2[2]], '-', color=self.color)
