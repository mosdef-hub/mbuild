__author__ = 'sallai'
import numpy as np
from mbuild.atom import *


class Angle(object):
    @classmethod
    def create(cls, atom1, atom2, atom3, angleType='undefined', color='black'):
        b = Angle()
        b.angleType = angleType
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

    def __hash__(self):
        # atom1 and atom3 are interchangeable
        return ( self.atom1.__hash__() * self.atom3.__hash__()) ^ self.atom1.__hash__() ^ self.atom3.__hash__() ^ self.atom2.__hash__()

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def plot(self, ax):
        epsilon = .4
        offset = np.array([round(random()*10+1)*.05, 0, 0])

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
