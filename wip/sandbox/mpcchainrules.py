from mbuild.prototype import Prototype

__author__ = 'sallai'

from mbuild.rules import RuleEngine
from mbuild.file_formats.xyz import *
from mbuild.angle import *
from mbuild.dihedral import *


class MpcChainRules(RuleEngine):

    def __init__(self, compound):
        super(MpcChainRules, self).__init__(compound)

    def execute(self):

        self.add_bond("CB", "HB", .9, 1.4, "c-h")
        self.add_bond("CB", "CB", 1.0, 1.6, "c-c")
        self.add_angle("HB", "CB", "CB", "h-c-c")

        # backbone
        self.add_bond("CB", "HB", .9, 1.4, "c-h")
        self.add_bond("CB", "CB", 1.0, 1.6, "c-c")

        self.add_angle("HB", "CB", "CB", "h-c-c")
        self.add_angle("HB", "CB", "HB", "h-c-h")
        self.add_angle("CB", "CB", "CB", "c-c-c")

        self.add_dihedral("HB", "CB", "CB", "CB", "h-c-c-c")
        self.add_dihedral("CB", "CB", "CB", "CB", "c-c-c-c")

        #link
        self.add_bond("CB", "C", 1.0, 2.2, "c-c")

        # chain
        self.add_bond("C", "H", .9, 1.4, "c-h")
        self.add_bond("C2", "H", .9, 1.4, "c2-h")
        self.add_bond("C", "C2", 1.0, 1.4, "c-c2")
        self.add_bond("C", "O", 1.2, 1.6, "c-o")
        self.add_bond("C", "C", 1.0, 1.8, "c-c")
        self.add_bond("N", "C", 1.2, 1.6, "n-c")
        self.add_bond("P", "O", 1.2, 1.9, "p-o")

        self.add_angle("C", "C", "O", "c-c-o")
        self.add_angle("H", "C", "C", "h-c-c")
        self.add_angle("C", "C", "C", "c-c-c")
        self.add_angle("H", "C2", "C", "h-c2-c")
        self.add_angle("C2", "C", "C", "c2-c-c")
        self.add_angle("H", "C", "H", "h-c-h")
        self.add_angle("C", "O", "C", "c-o-c")
        self.add_angle("C", "N", "C", "c-n-c")
        self.add_angle("O", "P", "O", "o-p-o")
        self.add_angle("C", "O", "P", "c-o-p")
        self.add_angle("O", "C", "O", "o-c-o")

        # # self.add_dihedral(H, C, C, C, "hxcxcxc", color=(0, 0, 0))
        self.add_dihedral("H", "C2", "C", "C", "h-c2-c-c")

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

    Prototype("c-h", color=(.8, .8, .8))
    Prototype("c-c", color=(0, .8, 0))

    Prototype("h-c-c", color=(1, 0.5, 0))

    # backbone
    Prototype("c-h", color=(.8, .8, .8))
    Prototype("c-c", color=(0, .8, 0))

    Prototype("h-c-c", color=(1, 0.5, 0))
    Prototype("h-c-h", color=(1, 0.5, 0))
    Prototype("c-c-c", color=(1, 0, 0))

    Prototype("h-c-c-c", color=(0, 0, 0))
    Prototype("c-c-c-c", color=(0.2, 0, 0))

    #link
    Prototype("c-c", color=(1, 0, 0))

    # chain
    Prototype("c-h", color=(1, 1, 1))
    Prototype("c2-h", color=(.8, .8, .8))
    Prototype("c-c2", color=(0, .8, 0))
    Prototype("c-o", color=(0, 1, 0))
    Prototype("c-c", color=(0, 1, 0))
    Prototype("n-c", color=(0, 0, 1))
    Prototype("p-o", color=(1, 0, 1))

    Prototype("c-c-o", color=(1, 0.5, 0))
    Prototype("h-c-c", color=(1, 0.5, 0))
    Prototype("c-c-c", color=(1, 0, 0))
    Prototype("h-c2-c", color=(1, 0.5, 0))
    Prototype("c2-c-c", color=(.8, 0, 0))
    Prototype("h-c-h", color=(1, 0.5, 0))
    Prototype("c-o-c", color=(1, 0, 0.5))
    Prototype("c-n-c", color=(0, 0, 0.5))
    Prototype("o-p-o", color=(.8, 0, 0.8))
    Prototype("c-o-p", color=(.8, 0.8, 0))
    Prototype("o-c-o", color=(.3, 0.8, 0))

    Prototype("h-c2-c-c", color=(0, 0, 0))


    Plot(mm).show()
