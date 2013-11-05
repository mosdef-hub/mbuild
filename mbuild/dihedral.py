__author__ = 'sallai'
import numpy as np
import pdb
from mbuild.atom import *


class Dihedral(object):
    @classmethod
    def create(cls, atom1, atom2, atom3, atom4, dihedralType='undefined', color='black'):
        b = dihedral()
        b.dihedralType = dihedralType
        b.color = color
        b.atom2 = atom2
        if atom1.__hash__() < atom2.__hash__():
            b.atom1 = atom1
            b.atom3 = atom3
        else:
            b.atom1 = atom3
            b.atom3 = atom1
        return b

    @classmethod
    def createFromAngles(cls, angle1, angle2, angle3, **kwargs):
        if (angle1.atom1, angle1.atom2) == (angle2.atom2, angle2.atom3):
            return Dihedral.create(angle1.atom3, angle1.atom2, angle2.atom2, angle2.atom1, **kwargs)
        elif (angle1.atom1, angle1.atom2) == (angle2.atom1, angle2.atom2):
            return Dihedral.create(angle1.atom3, angle1.atom2, angle2.atom2, angle2.atom3, **kwargs)
        elif (angle1.atom2, angle1.atom3) == (angle2.atom1, angle2.atom2):
            return Dihedral.create(angle1.atom1, angle1.atom2, angle2.atom2, angle2.atom3, **kwargs)
        elif (angle1.atom2, angle1.atom3) == (angle2.atom2, angle2.atom3):
            return Dihedral.create(angle1.atom1, angle1.atom2, angle2.atom2, angle2.atom1, **kwargs)
        else:
            return None

    """
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
        try:
            return math.acos(dot(a, b) / (linalg.norm(a) * linalg.norm(b)))
        except:
            pdb.set_trace()

    @staticmethod
    def computeInDegrees(atom1, atom2, atom3):
        return Angle.computeInRadians(atom1, atom2, atom3) * 180 / math.pi


    def inRadians(self):
        return Angle.computeInRadians(self.atom1, self.atom2, self.atom3)

    def inDegrees(self):
        return self.inRadians() * 180 / math.pi
    """

    def __hash__(self):
        # atom1 and atom3 are interchangeable
        return self.atom1.__hash__() ^ self.atom2.__hash__() ^ self.atom3.__hash__() ^ self.atom4.__hash__()


    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def plot(self, ax):
        epsilon = .4
        offset = np.array([np.round(np.random.random()*10+1)*.05, 0, 0])

#        offset = np.array([round(random()*10+1)*.05, 0, 0])

        pos1 = np.array(self.atom1.pos)
        pos2 = np.array(self.atom2.pos)
        pos3 = np.array(self.atom3.pos)
        v21 = pos1 - pos2 # vector from atom2 to atom1
        d21 = np.linalg.norm(v21) # atom1-atom2 distance
        v23 = pos3 - pos2 # vector from atom2 to atom3
        d23 = np.linalg.norm(v23) # atom3-atom2 distance

        p2 = pos2 + offset
        p1 = pos2 + v21*epsilon + offset
        p3 = pos2 + v23*epsilon + offset

        ax.plot([p1[0], p2[0], p3[0]],
                [p1[1], p2[1], p3[1]],
                [p1[2], p2[2], p3[2]], '-', color=self.color)

    def hasTypes(self, atomType1, atomType2, atomType3, atomType4):
        abcd = Dihedral.orderDihedral(self, atomType1, atomType2, atomType3, atomType4)
        if abcd:
            return True
        else:
            return False

    @classmethod
    def orderDihedral(cls, dihedral, type_A, type_B, type_C, type_D):
        abcd = Dihedral()
        abcd.dihedralType = dihedral.dihedralType
        abcd.color = dihedral.color
        if (isinstance(dihedral.atom1, type_A) and
            isinstance(dihedral.atom2, type_B) and
            isinstance(dihedral.atom3, type_C) and
            isinstance(dihedral.atom4, type_D)):
            abcd.atom1 = dihedral.atom1
            abcd.atom2 = dihedral.atom2
            abcd.atom3 = dihedral.atom3
            abcd.atom4 = dihedral.atom4
            return abcd
        elif (isinstance(dihedral.atom1, type_D) and
            isinstance(dihedral.atom2, type_C) and
            isinstance(dihedral.atom3, type_B) and
            isinstance(dihedral.atom4, type_A)):
            abcd.atom1 = dihedral.atom4
            abcd.atom2 = dihedral.atom3
            abcd.atom3 = dihedral.atom2
            abcd.atom4 = dihedral.atom1
            return abcd
