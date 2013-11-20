from __builtin__ import classmethod
from itertools import *
from mbuild.moleculemodel import MoleculeModel

from mbuild.xyz import *
from mbuild.bond import *
from mbuild.angle import *
from mbuild.dihedral import *

__author__ = 'sallai'


class RuleEngine(object):
    dmin = 0.9
    dmax = 1.1

    @classmethod
    def create(cls, model):
        re = RuleEngine()
        re.model = model
        return re


    def execute(self):
        self.add_bond(C, H, .9, 1.4, "c-h", (1, 1, 1))
        self.add_bond(C2, H, .9, 1.4, "c2-h", (.8, .8, .8))
        self.add_bond(C, C2, 1.0, 1.4, "c-c2", (0, .8, 0))
        self.add_bond(C, O, 1.2, 1.6, "c-o", (0, 1, 0))
        self.add_bond(C, C, 1.0, 1.8, "c-c", (0, 1, 0))
        self.add_bond(N, C, 1.2, 1.6, "n-c", (0, 0, 1))
        self.add_bond(P, O, 1.2, 1.9, "p-o", (1, 0, 1))

        self.add_angle(C, C, O, "c-c-o", color=(1, 0.5, 0))
        self.add_angle(H, C, C, "h-c-c", color=(1, 0.5, 0))
        self.add_angle(C, C, C, "c-c-c", color=(1, 0, 0))
        self.add_angle(H, C2, C, "h-c2-c", color=(1, 0.5, 0))
        self.add_angle(C2, C, C, "c2-c-c", color=(.8, 0, 0))
        self.add_angle(H, C, H, "h-c-h", color=(1, 0.5, 0))
        self.add_angle(C, O, C, "c-o-c", color=(1, 0, 0.5))
        self.add_angle(C, N, C, "c-n-c", color=(0, 0, 0.5))
        self.add_angle(O, P, O, "o-p-o", color=(.8, 0, 0.8))
        self.add_angle(C, O, P, "c-o-p", color=(.8, 0.8, 0))
        self.add_angle(O, C, O, "o-c-o", color=(.3, 0.8, 0))

        # self.add_dihedral(H, C, C, C, "hxcxcxc", color=(0, 0, 0))
        self.add_dihedral(H, C2, C, C, "h-c2-c-c", color=(0, 0, 0))

    def add_bond(self, type_A, type_B, dmin, dmax, kind, color=(1,1,1)):
        "Ai-Bj distance is in [dmin, dmax] => add bond A1xB(Ai,Bj) (symmetric)"
        for a1 in self.model.getAtomsByKind(type_A):
            nearest = self.model.getAtomsInRange(a1.pos, 2)
            for b1 in nearest:
                if isinstance(b1, type_B) and (dmin <= b1.distance(a1) <= dmax):
                    self.model.add(Bond.create(a1, b1, kind=kind, color=color))

    def add_angle(self, type_A, type_B, type_C, kind, thmin=-Inf, thmax=Inf, color=(1,1,1)):
        """
        """
        for ab1 in self.model.getBondsByAtomKind(type_A, type_B):
            ab = ab1.cloneWithOrder(type_A, type_B)
            nearest = self.model.getBondsInRange(ab.com(), 3)
            for bc1 in nearest:
                if ab1 == bc1:
                    continue
                if not bc1.hasAtomKinds(type_B, type_C):
                    continue

                bc = bc1.cloneWithOrder(type_B, type_C)

                temp_ang = Angle.createFromBonds(ab, bc, kind=kind, color=color)
                if temp_ang:
                    if isinstance(temp_ang.atom2, type_B) and (thmin <= temp_ang.inDegrees() <= thmax):
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

    def c1xc1xc(self):
        for b1, b2 in ifilter(lambda (bond1, bond2): bond1.hasCommonAtomsWith(bond2),
                              combinations(
                                      ifilter(lambda bond: bond.kind == 'c1xc',
                                              self.model.bonds
                                      ),
                                      2
                              )
        ):
            # create angle
            a = Angle.createFromBonds(b1, b2)
            angle_deg = a.inDegrees()
            if 106 <= angle_deg <= 113: # around 109.5
                a.kind = 'c1xc1xc108'
                a.color = 'red'
                self.model.add(a)
            elif 118 <= angle_deg <= 122: # around 120
                a.kind = 'c1xc1xc120'
                a.color = 'blue'
                self.model.add(a)
            else:
                print "bad angle:" + str(angle_deg)


    def c1xc1xc_simple(self):
        # Ci-Cj-Ck and Ci != Ck => add angle C1xC1xC(Ci,Cj,Ck) (symmetric)
        amin = 0
        amax = 1

        # for all c
        for ci in self.model.atoms:
            if isinstance(ci, C):
                # for all c1xc bonds that ci is part of
                for ci1xcj in ci.bonds:
                    if ci1xcj.kind == 'c1xc':
                        cj = ci1xcj.atom1 if ci1xcj.atom1 != ci else ci1xcj.atom2
                        # for all other c1xc bonds of ci
                        for ci1xck in ci.bonds:
                            if ci1xck.kind == 'c1xc' and ci1xck.atom1 != cj and ci1xck.atom2 != cj:
                                ck = ci1xck.atom1 if ci1xck.atom1 != ci else ci1xck.atom2
                                # create angle
                                angle_deg = Angle.computeInDegrees(cj, ci, ck)
                                if 106 <= angle_deg <= 113:
                                    a = Angle.create(cj, ci, ck, 'c1xc1xc108')
                                    a.color = 'red'
                                    self.model.add(a)
                                elif 118 <= angle_deg <= 122:
                                    a = Angle.create(cj, ci, ck, 'c1xc1xc120')
                                    a.color = 'blue'
                                    self.model.add(a)
                                else:
                                    print "bad angle:" + str(angle_deg)


                                    # "Ci is bonded to Cj and Ci is bonded to Hk => add angle CC1xH(Cj,Ci,Hk) (symmetric)"


                                    # Hi-Cj-Ck-Hl => add dihedral HCC1xH(Hi,Cj,Ck,Hl) (symmetric)


if __name__ == "__main__":
    # m = Xyz.create("c60.xyz")
    m = Xyz.create("../mpc.xyz")
    mm = MoleculeModel.create()
    mm.add([atom for label, atom in m.atoms()])


    r = RuleEngine.create(mm)

    r.execute()
    print 'n_bonds: ' + str(len(r.model.bonds))
    print 'n_angles: ' + str(len(r.model.angles))
    print 'n_dihedrals: ' + str(len(r.model.dihedrals))

    # for a in r.model.angles:
    #     print a.inDegrees()

    #print(mm.getAtomsInRange(mm.atoms.pop().pos,2))
    print "Missing angle kinds: " + str(mm.findMissingAngleKinds())
    r.model.plot(angles=True)


    # m.plot(labels=False)


