from mbuild.plot import Plot

__author__ = 'sallai'

from __builtin__ import classmethod
from itertools import *

from mbuild.moleculemodel import MoleculeModel
from mbuild.rules import RuleEngine
from mbuild.xyz import *
from mbuild.bond import *
from mbuild.angle import *
from mbuild.dihedral import *
from mbuild.mdio import write_lammps_data



class MpcChainRules(RuleEngine):

    def __init__(self, model):
        super(MpcChainRules, self).__init__(model)

    def execute(self):

        self.add_bond("CB", "HB", .9, 1.4, "c-h", color=(.8, .8, .8))
        self.add_bond("CB", "CB", 1.0, 1.6, "c-c", color=(0, .8, 0))
        self.add_angle("HB", "CB", "CB", "h-c-c", color=(1, 0.5, 0))

        # backbone
        self.add_bond("CB", "HB", .9, 1.4, "c-h", (.8, .8, .8))
        self.add_bond("CB", "CB", 1.0, 1.6, "c-c", (0, .8, 0))

        self.add_angle("HB", "CB", "CB", "h-c-c", color=(1, 0.5, 0))
        self.add_angle("HB", "CB", "HB", "h-c-h", color=(1, 0.5, 0))
        self.add_angle("CB", "CB", "CB", "c-c-c", color=(1, 0, 0))

        self.add_dihedral("HB", "CB", "CB", "CB", "h-c-c-c", color=(0, 0, 0))
        self.add_dihedral("CB", "CB", "CB", "CB", "c-c-c-c", color=(0.2, 0, 0))

        #link
        self.add_bond("CB", "C", 1.0, 2.2, "c-c", (1, 0, 0))

        # chain
        self.add_bond("C", "H", .9, 1.4, "c-h", (1, 1, 1))
        self.add_bond("C2", "H", .9, 1.4, "c2-h", (.8, .8, .8))
        self.add_bond("C", "C2", 1.0, 1.4, "c-c2", (0, .8, 0))
        self.add_bond("C", "O", 1.2, 1.6, "c-o", (0, 1, 0))
        self.add_bond("C", "C", 1.0, 1.8, "c-c", (0, 1, 0))
        self.add_bond("N", "C", 1.2, 1.6, "n-c", (0, 0, 1))
        self.add_bond("P", "O", 1.2, 1.9, "p-o", (1, 0, 1))

        self.add_angle("C", "C", "O", "c-c-o", color=(1, 0.5, 0))
        self.add_angle("H", "C", "C", "h-c-c", color=(1, 0.5, 0))
        self.add_angle("C", "C", "C", "c-c-c", color=(1, 0, 0))
        self.add_angle("H", "C2", "C", "h-c2-c", color=(1, 0.5, 0))
        self.add_angle("C2", "C", "C", "c2-c-c", color=(.8, 0, 0))
        self.add_angle("H", "C", "H", "h-c-h", color=(1, 0.5, 0))
        self.add_angle("C", "O", "C", "c-o-c", color=(1, 0, 0.5))
        self.add_angle("C", "N", "C", "c-n-c", color=(0, 0, 0.5))
        self.add_angle("O", "P", "O", "o-p-o", color=(.8, 0, 0.8))
        self.add_angle("C", "O", "P", "c-o-p", color=(.8, 0.8, 0))
        self.add_angle("O", "C", "O", "o-c-o", color=(.3, 0.8, 0))

        # # self.add_dihedral(H, C, C, C, "hxcxcxc", color=(0, 0, 0))
        self.add_dihedral("H", "C2", "C", "C", "h-c2-c-c", color=(0, 0, 0))

if __name__ == "__main__":
    # m = Xyz.create("c60.xyz")
    mm = Xyz("mpcchain.xyz")

    MpcChainRules(mm).execute()

    print 'n_bonds: ' + str(len(mm.bonds))
    print 'n_angles: ' + str(len(mm.angles))
    print 'n_dihedrals: ' + str(len(mm.dihedrals))


    for bond in mm.bonds:
        print str(bond) +":"+ bond.kind + " "+ bond.atom1.kind + " " + bond.atom2.kind + " " + str(bond.atom1) + " " + str(bond.atom2)

    for bond in mm.getBondsByAtomKind("CB","HB"):
        print bond

    #print "Missing angle kinds: " + str(mm.findMissingAngleKinds())
    #r.model.plot(angles=True)

    # write_lammps_data(mm, 'data.mpcchain', 'pMPC')
    # pdb.set_trace()

    Plot(mm).show()

    print "hi"