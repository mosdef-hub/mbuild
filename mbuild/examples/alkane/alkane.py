__author__ = 'sallai'

from mbuild.compound import Compound
from mbuild.coordinate_transform import *
from mbuild.examples.alkane.ch2 import Ch2
from mbuild.examples.ethane.methyl import Methyl
from mbuild.polymer import Polymer


class Alkane(Compound):
    """ """
    def __init__(self, n=3):
        if n < 2:
            raise Exception('n must be 1 or more')
        Compound.__init__(self)

        self.add(Methyl(), "methyl1")
        self.add(Polymer(Ch2(), n=n-2, port_labels=("up", "down")), "chain")
        self.add(Methyl(), "methyl2")

        equivalence_transform(self.chain, self.chain.up, self.methyl1.down)
        equivalence_transform(self.methyl2, self.methyl2.up, self.chain.down)

if __name__ == "__main__":
    n = 3
    alkane = Alkane(n=n)

    #m.save("{}-alkane.pdb".format(n))

    # mol = m.to_molecule()
    # # mol.localopt(forcefield="gaff", steps=2000)
    # mol.localopt(forcefield="gaff", steps=2000)
    # print mol
    # m.update_from_molecule(mol)


    from mbuild.tools import add_angle

    alkane = alkane.to_trajectory()
    alkane.top.load_ff_bonds()

    #alkane.top.enumerate_ff_angles()
    #alkane.top.enumerate_ff_dihedrals()
    alkane.top.enumerate_ff_angles_and_dihedrals()

    print len(alkane.top._ff_bonds)
    print len(alkane.top._ff_angles)
    print len(alkane.top._ff_dihedrals)
    import pdb
    pdb.set_trace()

    from mbuild.plot import Plot
    Plot(alkane).show()
