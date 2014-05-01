from warnings import warn

__author__ = 'sallai'
# import numpy as np
# import pdb
from atom import *
import copy

class Angle(object):

    def __init__(self, atom1, atom2, atom3, kind='undefined'):
        assert(isinstance(kind, basestring))
        assert(not atom1 == atom2)
        assert(not atom2 == atom3)
        assert(not atom3 == atom1)
        self.kind = kind
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

    @classmethod
    def createFromBonds(cls, bond1, bond2, **kwargs):
        if bond1.atom1 == bond2.atom1:
            return Angle(bond1.atom2, bond1.atom1, bond2.atom2, **kwargs)
        if bond1.atom1 == bond2.atom2:
            return Angle(bond1.atom2, bond1.atom1, bond2.atom1, **kwargs)
        elif bond1.atom2 == bond2.atom1:
            return Angle(bond1.atom1, bond1.atom2, bond2.atom2, **kwargs)
        elif bond1.atom2 == bond2.atom2:
            return Angle(bond1.atom1, bond1.atom2, bond2.atom1, **kwargs)
        else:
            return None

    def clone(self):
        return copy.copy(self)

    def cloneImage(self):
        new_obj = self.clone()
        new_obj.atom1 = self.atom3
        new_obj.atom2 = self.atom2
        new_obj.atom3 = self.atom1
        return new_obj

    def cloneWithOrder(self, type_A, type_B, type_C):
        assert(isinstance(type_A, basestring))
        assert(isinstance(type_B, basestring))
        assert(isinstance(type_C, basestring))
        if (self.atom1.kind==type_A) and (self.atom2.kind==type_B)  and (self.atom3.kind==type_C):
            return self.clone()
        elif (self.atom1.kind==type_C) and (self.atom2.kind==type_B) and (self.atom3.kind==type_A):
            return self.cloneImage()

        warn("cannot clone angle " + str(self) + " with order " + str(type_A) + "," + str(type_B)+ "," + str(type_C))


    @staticmethod
    def computeInRadians(atom1, atom2, atom3):
        """
        Compute the angle at atom2
        :param atom1:
        :param atom2:
        :param atom3:
        :return: the angle in radians
        """
        assert isinstance(atom1, Atom)
        assert isinstance(atom2, Atom)
        assert isinstance(atom3, Atom)
        a = array(atom1.pos) - array(atom2.pos)
        b = array(atom3.pos) - array(atom2.pos)
        return math.acos(dot(a, b) / (linalg.norm(a) * linalg.norm(b)))

    @staticmethod
    def computeInDegrees(atom1, atom2, atom3):
        return Angle.computeInRadians(atom1, atom2, atom3) * 180 / math.pi

    def inRadians(self):
        return Angle.computeInRadians(self.atom1, self.atom2, self.atom3)

    def inDegrees(self):
        return self.inRadians() * 180 / math.pi

    def hasAtomKinds(self, atomType1, atomType2, atomType3):
        """
        if isinstance(atomType1, type):
            atomType1 = atomType1.kind
        if isinstance(atomType2, type):
            atomType2 = atomType2.kind
        if isinstance(atomType3, type):
            atomType3 = atomType3.kind
        """
        return ((self.atom1.kind == atomType1 and self.atom2.kind == atomType2 and self.atom3.kind == atomType3) or
                (self.atom1.kind == atomType3 and self.atom2.kind == atomType2 and self.atom3.kind == atomType1))


    def __repr__(self):
        return "Angle"+str(id(self))+"("+ str(self.atom1) + ", " + str(self.atom2) + ", " + str(self.atom3) + ", kind=" + self.kind+")"

    def __hash__(self):
        return hash((self.atom1, self.atom2, self.atom3))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
