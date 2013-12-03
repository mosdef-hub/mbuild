from __builtin__ import classmethod
from itertools import *
from mbuild.moleculemodel import MoleculeModel
from mbuild.rules import RuleEngine

from mbuild.xyz import *
from mbuild.bond import *
from mbuild.angle import *
from mbuild.dihedral import *

__author__ = 'sallai'


class MpcChainRules(RuleEngine):

    @classmethod
    def create(cls, model):
        re = super(MpcChainRules, cls).create(model)
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


if __name__ == "__main__":
    # m = Xyz.create("c60.xyz")
    m = Xyz.create("../mpcchain.xyz")
    mm = MoleculeModel.create()
    mm.add([atom for label, atom in m.atoms()])


    r = MpcChainRules.create(mm)

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


