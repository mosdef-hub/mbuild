from __builtin__ import classmethod
from itertools import *
from mbuild.moleculemodel import MoleculeModel

from mbuild.xyz import *
from mbuild.bond import *
from mbuild.angle import *
from mbuild.dihedral import *

__author__ = 'sallai'


class RuleEngine(object):

    def __init__(self, model):
        self.model = model

    def execute(self):
        raise Exception("RuleEngine must be subclassed, with the execute method implemented in the subclass")

    def add_bond(self, type_A, type_B, dmin, dmax, kind, color=(1.0,1.0,1.0)):
        "Ai-Bj distance is in [dmin, dmax] => add bond A1xB(Ai,Bj) (symmetric)"
        for a1 in self.model.getAtomsByKind(type_A):
            nearest = self.model.getAtomsInRange(a1.pos, dmax)
            for b1 in nearest:
                if (b1.kind==type_B) and (dmin <= b1.distance(a1) <= dmax):
                    self.model.add(Bond(a1, b1, kind=kind, color=color))

    def add_angle(self, type_A, type_B, type_C, kind, thmin=-Inf, thmax=Inf, color=(1,1,1)):
        """
        """

        for ab1 in self.model.getBondsByAtomKind(type_A, type_B):
            ab = ab1.cloneWithOrder(type_A, type_B)
            nearest = self.model.getBondsInRange(ab.com(), 10)
            for bc1 in nearest:
                if ab1 == bc1:
                    continue
                if not bc1.hasAtomKinds(type_B, type_C):
                    continue

                bc = bc1.cloneWithOrder(type_B, type_C)

                temp_ang = Angle.createFromBonds(ab, bc, kind=kind, color=color)
                if temp_ang:
                    if (temp_ang.atom2.kind==type_B) and (thmin <= temp_ang.inDegrees() <= thmax):
                        self.model.add(temp_ang)

    def add_dihedral(self, type_A, type_B, type_C, type_D, dihedralKind, color=(1,1,1)):
        """
        """

        for abc1 in self.model.getAnglesByAtomKind(type_A, type_B, type_C):
            abc = abc1.cloneWithOrder(type_A, type_B, type_C)

            nearest = self.model.getAnglesInRange(abc.atom2.pos, 5)
            for bcd1 in nearest:
                if abc1 == bcd1:
                    continue

                if not bcd1.hasAtomKinds(type_B, type_C, type_D):
                    continue

                bcd = bcd1.cloneWithOrder(type_B, type_C, type_D)

                temp_dhdr = Dihedral.createFromAngles(abc, bcd, kind=dihedralKind, color=color)
                if temp_dhdr:
                    self.model.add(temp_dhdr)
