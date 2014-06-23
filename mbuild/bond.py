from copy import copy, deepcopy
from warnings import warn
# from mbuild.atom import Atom
# import numpy as np

__author__ = 'sallai'
class Bond(object):

    __slots__ = ['atom1','atom2', 'kind']
    def __init__(self, atom1, atom2, kind='undefined'):
        assert(isinstance(kind, basestring))
        assert(not atom1 == atom2)
        self.atom1 = atom1
        self.atom2 = atom2
        self.kind = kind

    def clone(self):
        return copy(self)

    def cloneImage(self):
        b = self.clone()
        b.atom1 = self.atom2
        b.atom2 = self.atom1
        return b

    def cloneWithOrder(self, type_A, type_B):
        ab = self.clone()
        if (self.atom1.kind==type_A) and (self.atom2.kind==type_B):
            ab.atom1 = self.atom1
            ab.atom2 = self.atom2
            return ab
        elif (self.atom1.kind==type_B) and (self.atom2.kind==type_A):
            ab.atom1 = self.atom2
            ab.atom2 = self.atom1
            return ab
        warn("cannot clone bond " + str(self) + " with order " + str(type_A) + "," + str(type_B))

    @classmethod
    def orderBond(cls, bond, type_A, type_B):
        cls.cloneWithOrder(bond, type_A, type_B)

    def com(self):
        return [sum(y) / 2 for y in zip(self.atom1.pos, self.atom2.pos)]

    def hasCommonAtomsWith(self, other):
        assert(isinstance(other, Bond))
        if self.atom1 == other.atom1 or self.atom1 == other.atom2 or self.atom2 == other.atom1 or self.atom2 == other.atom2:
            return True
        else:
            return False

    def hasAtomKinds(self, atomType1, atomType2):
        return ((self.atom1.kind == atomType1 and self.atom2.kind == atomType2) or
                (self.atom1.kind == atomType2 and self.atom2.kind == atomType1))

    def __hash__(self):
        return hash((self.atom1, self.atom2))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __repr__(self):
        return "Bond"+str(id(self))+"("+str(self.atom1)+","+str(self.atom2)+", kind="+self.kind+")"


    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        memo[id(self)] = newone

        newone.kind = deepcopy(self.kind, memo)
        newone.atom1 = deepcopy(self.atom1, memo)
        newone.atom2 = deepcopy(self.atom2, memo)

        # for k, v in self.__dict__.items():
        #     setattr(newone, k, deepcopy(v, memo))

        return newone

