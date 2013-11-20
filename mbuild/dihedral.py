__author__ = 'sallai'
import numpy as np
import pdb
from mbuild.atom import *
import copy

class Dihedral(object):
    @classmethod
    def create(cls, atom1, atom2, atom3, atom4, kind='undefined', color='black'):
        b = Dihedral()
        b.kind = kind
        b.color = color
        b.atom1 = atom1
        b.atom2 = atom2
        b.atom3 = atom3
        b.atom4 = atom4
        return b

    @classmethod
    def createFromAngles(cls, angle1, angle2, **kwargs):
        if (angle1.atom1, angle1.atom2) == (angle2.atom2, angle2.atom3):
            return Dihedral.create(angle1.atom3, angle1.atom2, angle2.atom2, angle2.atom1, **kwargs)
        elif (angle1.atom1, angle1.atom2) == (angle2.atom2, angle2.atom1):
            return Dihedral.create(angle1.atom3, angle1.atom2, angle2.atom2, angle2.atom3, **kwargs)
        elif (angle1.atom2, angle1.atom3) == (angle2.atom1, angle2.atom2):
            return Dihedral.create(angle1.atom1, angle1.atom2, angle2.atom2, angle2.atom3, **kwargs)
        elif (angle1.atom2, angle1.atom3) == (angle2.atom3, angle2.atom2):
            return Dihedral.create(angle1.atom1, angle1.atom2, angle2.atom2, angle2.atom1, **kwargs)
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

        warn ("cannot clone dihedral " + str(self) + " with order " + str(type_A) + "," + str(type_B)+ "," + str(type_C)+ "," + str(type_D))

    def __repr__( self ):
        return "Dihedral[atom1=" + str(self.atom1) + " atom2=" + str(self.atom2) + " atom3=" + str(self.atom3) + " atom4=" + str(self.atom4) + "]"

    def __hash__(self):
        # return hash((self.kind, self.atom1, self.atom2, self.atom3, self.atom4))
        return hash((self.atom1, self.atom2, self.atom3, self.atom4))


    # def __hash__(self):
    #     # atom1 and atom3 are interchangeable
    #     return self.atom1.__hash__() ^ self.atom2.__hash__() ^ self.atom3.__hash__() ^ self.atom4.__hash__()
    #
    #
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

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

if __name__ == "__main__":
    atom1 = C(pos=(0,-1,0))
    atom2 = C(pos=(0,0,0))
    atom3 = C(pos=(1,0,0))
    atom4 = C(pos=(1,1,1))
    d = Dihedral.create(atom1, atom2, atom3, atom4)
    print "Dihedral angle should be 135 degrees, we compute it to be: " + str(d.inDegrees())

