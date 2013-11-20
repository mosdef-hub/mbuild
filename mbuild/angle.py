from warnings import warn

__author__ = 'sallai'
import numpy as np
import pdb
from mbuild.atom import *
import copy

class Angle(object):
    @classmethod
    def create(cls, atom1, atom2, atom3, kind='undefined', color='black'):
        b = Angle()
        b.kind = kind
        b.color = color
        b.atom1 = atom1
        b.atom2 = atom2
        b.atom3 = atom3

        return b

    @classmethod
    def createFromBonds(cls, bond1, bond2, **kwargs):
        if bond1.atom1 == bond2.atom1:
            return Angle.create(bond1.atom2, bond1.atom1, bond2.atom2, **kwargs)
        if bond1.atom1 == bond2.atom2:
            return Angle.create(bond1.atom2, bond1.atom1, bond2.atom1, **kwargs)
        elif bond1.atom2 == bond2.atom1:
            return Angle.create(bond1.atom1, bond1.atom2, bond2.atom2, **kwargs)
        elif bond1.atom2 == bond2.atom2:
            return Angle.create(bond1.atom1, bond1.atom2, bond2.atom1, **kwargs)
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
        if isinstance(self.atom1, type_A) and isinstance(self.atom2, type_B)  and isinstance(self.atom3, type_C):
            return self.clone()
        elif isinstance(self.atom1, type_C) and isinstance(self.atom2, type_B) and isinstance(self.atom3, type_A):
            return self.cloneImage()

        warn ("cannot clone angle " + str(self) + " with order " + str(type_A) + "," + str(type_B)+ "," + str(type_C))


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
        if isinstance(atomType1, type):
            atomType1 = atomType1.kind
        if isinstance(atomType2, type):
            atomType2 = atomType2.kind
        if isinstance(atomType3, type):
            atomType3 = atomType3.kind

        return (self.atom1.kind == atomType1 and self.atom2.kind == atomType2 and self.atom3.kind == atomType3) or (self.atom1.kind == atomType3 and self.atom2.kind == atomType2 and self.atom3.kind == atomType1)


    def __repr__( self ):
        return "Angle[atom1=" + str(self.atom1) + " atom2=" + str(self.atom2) + " atom3=" + str(self.atom3) + "]"


    # def __hash__(self):
    #     # atom1 and atom3 are interchangeable
    #     return ( self.atom1.__hash__() * self.atom3.__hash__()) ^ self.atom1.__hash__() ^ self.atom3.__hash__() ^ self.atom2.__hash__()
    #
    # def __eq__(self, other):
    #     return self.__hash__() == other.__hash__()

#     def plot(self, ax):
#         epsilon = .4
#         offset = np.array([np.round(np.random.random()*10+1)*.05, 0, 0])
#
# #        offset = np.array([round(random()*10+1)*.05, 0, 0])
#
#         pos1 = np.array(self.atom1.pos)
#         pos2 = np.array(self.atom2.pos)
#         pos3 = np.array(self.atom3.pos)
#         v21 = pos1 - pos2 # vector from atom2 to atom1
#         d21 = np.linalg.norm(v21) # atom1-atom2 distance
#         v23 = pos3 - pos2 # vector from atom2 to atom3
#         d23 = np.linalg.norm(v23) # atom3-atom2 distance
#
#         p2 = pos2 + offset
#         p1 = pos2 + v21*epsilon + offset
#         p3 = pos2 + v23*epsilon + offset
#
#         ax.plot([p1[0], p2[0], p3[0]],
#                 [p1[1], p2[1], p3[1]],
#                 [p1[2], p2[2], p3[2]], '-', color=self.color)

