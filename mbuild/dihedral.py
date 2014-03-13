__author__ = 'sallai'
import numpy as np
# import pdb
from atom import *
import copy

class Dihedral(object):

    def __init__(self, atom1, atom2, atom3, atom4, kind='undefined'):
        assert(isinstance(kind, basestring))
        assert(not atom1 == atom2)
        assert(not atom2 == atom3)
        assert(not atom3 == atom4)
        self.kind = kind
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    @classmethod
    def createFromAngles(cls, angle1, angle2, **kwargs):
        if (angle1.atom1, angle1.atom2) == (angle2.atom2, angle2.atom3):
            return Dihedral(angle1.atom3, angle1.atom2, angle2.atom2, angle2.atom1, **kwargs)
        elif (angle1.atom1, angle1.atom2) == (angle2.atom2, angle2.atom1):
            return Dihedral(angle1.atom3, angle1.atom2, angle2.atom2, angle2.atom3, **kwargs)
        elif (angle1.atom2, angle1.atom3) == (angle2.atom1, angle2.atom2):
            return Dihedral(angle1.atom1, angle1.atom2, angle2.atom2, angle2.atom3, **kwargs)
        elif (angle1.atom2, angle1.atom3) == (angle2.atom3, angle2.atom2):
            return Dihedral(angle1.atom1, angle1.atom2, angle2.atom2, angle2.atom1, **kwargs)
        else:
            return None

    def clone(self):
        return copy.copy(self)

    def cloneImage(self):
        new_obj = self.clone()
        new_obj.atom1 = self.atom4
        new_obj.atom2 = self.atom3
        new_obj.atom3 = self.atom2
        new_obj.atom4 = self.atom1
        return new_obj

    def cloneWithOrder(self, type_A, type_B, type_C, type_D):
        if isinstance(self.atom1, type_A) and isinstance(self.atom2, type_B) and isinstance(self.atom3, type_C) and isinstance(self.atom4, type_D):
            return self.clone()
        elif isinstance(self.atom1, type_D) and isinstance(self.atom2, type_C) and isinstance(self.atom3, type_B) and isinstance(self.atom4, type_A):
            return self.cloneImage()

        warn("cannot clone dihedral " + str(self) + " with order " + str(type_A) + "," + str(type_B)+ "," + str(type_C)+ "," + str(type_D))

    def hasAtomKinds(self, atomType1, atomType2, atomType3, atomType4):
        if isinstance(atomType1, type):
            atomType1 = atomType1.kind
        if isinstance(atomType2, type):
            atomType2 = atomType2.kind
        if isinstance(atomType3, type):
            atomType3 = atomType3.kind
        if isinstance(atomType4, type):
            atomType4 = atomType4.kind

        return (self.atom1.kind == atomType1 and self.atom2.kind == atomType2 and self.atom3.kind == atomType3 and self.atom4.kind == atomType4) or (self.atom1.kind == atomType4 and self.atom2.kind == atomType3 and self.atom3.kind == atomType2 and self.atom4.kind == atomType1)

    @staticmethod
    def computeNormal(atom1, atom2, atom3):
        # a is a vector from atom1 to atom 2
        a = np.array([atom2.pos[0]-atom1.pos[0],atom2.pos[1]-atom1.pos[1],atom2.pos[2]-atom1.pos[2]])
        # b is a vector from atom1 to atom 3
        b = np.array([atom3.pos[0]-atom1.pos[0],atom3.pos[1]-atom1.pos[1],atom3.pos[2]-atom1.pos[2]])
        # compute the cross product, which is perpendicular to both a and b, and thus, it's a normal vector to the plane
        print a
        print b
        return np.cross(a,b)

    @staticmethod
    def computeInRadians(atom1, atom2, atom3, atom4):
        assert isinstance(atom1, Atom)
        assert isinstance(atom2, Atom)
        assert isinstance(atom3, Atom)
        assert isinstance(atom4, Atom)

        # we need to compute the normal vectors of the plane defined by (atom1, atom2, atom3) and by (atom2, atom3, atom4)
        a = Dihedral.computeNormal(atom1, atom2, atom3)
        b = Dihedral.computeNormal(atom2, atom3, atom4)

        # then we compute and return the angle of those normals
        return math.acos(dot(a, b) / (linalg.norm(a) * linalg.norm(b)))

    @staticmethod
    def computeInDegrees(atom1, atom2, atom3, atom4):
        return Dihedral.computeInRadians(atom1, atom2, atom3, atom4) * 180 / math.pi

    def inRadians(self):
        return Dihedral.computeInRadians(self.atom1, self.atom2, self.atom3, self.atom4)

    def inDegrees(self):
        return self.inRadians() * 180 / math.pi

    # @classmethod
    # def orderDihedral(cls, dihedral, type_A, type_B, type_C, type_D):
    #     abcd = Dihedral()
    #     abcd.kind = dihedral.kind
    #     abcd.color = dihedral.color
    #     if (isinstance(dihedral.atom1, type_A) and
    #         isinstance(dihedral.atom2, type_B) and
    #         isinstance(dihedral.atom3, type_C) and
    #         isinstance(dihedral.atom4, type_D)):
    #         abcd.atom1 = dihedral.atom1
    #         abcd.atom2 = dihedral.atom2
    #         abcd.atom3 = dihedral.atom3
    #         abcd.atom4 = dihedral.atom4
    #         return abcd
    #     elif (isinstance(dihedral.atom1, type_D) and
    #         isinstance(dihedral.atom2, type_C) and
    #         isinstance(dihedral.atom3, type_B) and
    #         isinstance(dihedral.atom4, type_A)):
    #         abcd.atom1 = dihedral.atom4
    #         abcd.atom2 = dihedral.atom3
    #         abcd.atom3 = dihedral.atom2
    #         abcd.atom4 = dihedral.atom1
    #         return abcd

    def __repr__(self):
        return "Dihedral"+str(id(self))+"("+ str(self.atom1) + ", " + str(self.atom2) + ", " + str(self.atom3) + ", " + str(self.atom4) + ", kind=" + self.kind+")"

    def __hash__(self):
        return hash((self.atom1, self.atom2, self.atom3, self.atom4))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
