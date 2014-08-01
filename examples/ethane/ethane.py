from mbuild.bond import Bond

__author__ = 'sallai'

from methyl import Methyl
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.mol2file import load_mol2, write_mol2
from mbuild.coordinate_transform import *


class Ethane(Compound):
    """ """

    def __init__(self):
        super(Ethane, self).__init__(kind='Ethane')
        methyl1 = Methyl()
        methyl2 = Methyl()

        equivalence_transform(methyl2, methyl2.up, methyl1.down)

        self.add(methyl1, "m1")
        self.add(methyl2, "m2")

        # bond = Bond(methyl1.C_1, methyl2.C_1)
        bond = Bond(methyl1.up.C_1, methyl2.down.C_1)
        self.add(bond)


if __name__ == "__main__":
    ethane = Ethane()
    write_mol2(ethane, 'ethane.mol2')

    # import pdb
    # pdb.set_trace()

    from mbuild.plot import Plot
    Plot(ethane, verbose=True, atoms=True, bonds=True, angles=False, dihedrals=False).show()
