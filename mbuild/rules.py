__author__ = 'sallai'

import pdb

from xyz import *


class RuleEngine(object):

    def __init__(self, compound):
        self.compound = compound

    def execute(self):
        raise Exception("RuleEngine must be subclassed, with the execute method implemented in the subclass")

    def add_bond(self, type_A, type_B, dmin, dmax, kind):
        """Ai-Bj distance is in [dmin, dmax] => add bond A1xB(Ai,Bj) (symmetric)."""
        #for a1 in self.compound.getAtomsByBondType(type_A):
        for a1 in self.compound.getAtomListByKind(type_A):
            nearest = self.compound.getAtomsInRange(a1.pos, dmax, kind=type_B)
            for b1 in nearest:
                if (b1.kind==type_B) and (dmin <= self.compound.min_periodic_distance(b1.pos, a1.pos) <= dmax):
                    self.compound.add(Bond(a1, b1, kind=kind))
                    # print "Added bond of kind {0}".format(kind)

    def add_angle(self, type_A, type_B, type_C, kind):
        """
        """

        for ab1 in self.compound.getBondsByAtomKind(type_A, type_B):
            ab = ab1.cloneWithOrder(type_A, type_B)
            # nearest = self.compound.getBondsInRangeKdTree(ab.com(), 10)
            nearest = self.compound.getBondsInRange(ab)
            for bc1 in nearest:
                if ab1 == bc1:
                    continue
                if not bc1.hasAtomKinds(type_B, type_C):
                    continue

                bc = bc1.cloneWithOrder(type_B, type_C)

                temp_ang = Angle.createFromBonds(ab, bc, kind=kind)
                if (temp_ang) and (temp_ang.atom2.kind==type_B):
                    self.compound.add(temp_ang)

    def add_dihedral(self, type_A, type_B, type_C, type_D, dihedralKind):
        """
        """

        for abc1 in self.compound.getAnglesByAtomKind(type_A, type_B, type_C):
            abc = abc1.cloneWithOrder(type_A, type_B, type_C)

            nearest = self.compound.getAnglesInRange(abc.atom2.pos, 5)
            for bcd1 in nearest:
                if abc1 == bcd1:
                    continue

                if not bcd1.hasAtomKinds(type_B, type_C, type_D):
                    continue

                bcd = bcd1.cloneWithOrder(type_B, type_C, type_D)

                temp_dhdr = Dihedral.createFromAngles(abc, bcd, kind=dihedralKind)
                if temp_dhdr:
                    self.compound.add(temp_dhdr)
