from mbuild.plugins.system import System

__author__ = 'sallai'

from methyl import Methyl
from mbuild.compound import Compound
from mbuild.file_formats.mol2file import write_mol2
from mbuild.coordinate_transform import *


class Ethane(Compound):
    """ """

    def __init__(self):
        super(Ethane, self).__init__(kind='Ethane')
        self.add(Methyl(), "m1")
        self.add(Methyl(), "m2")

        equivalence_transform(self.m1, self.m1.up, self.m2.down)



if __name__ == "__main__":
    ethane = Ethane()
    write_mol2(ethane, 'ethane.mol2')

    # import pdb
    # pdb.set_trace()

    from mbuild.plot import Plot
    Plot(ethane, verbose=False, atoms=True, bonds=True, angles=False, dihedrals=False).show()
