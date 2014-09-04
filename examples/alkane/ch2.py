__author__ = 'sallai'
import os,sys
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.file_formats.mol2file import load_mol2
from mbuild.coordinate_transform import *

class Ch2(Compound):
    """ """
    def __init__(self):
        Compound.__init__(self)

        # Look for data file in same directory as this python module.
        self.append_from_file('ch2.pdb', relative_to_module=__name__)

        self.add(Port(anchor=self.C[0]), 'up')
        translate(self.up, np.array([0,0.7,0]))

        self.add(Port(anchor=self.C[0]), 'down')
        translate(self.down, np.array([0,-0.7,0]))

if __name__ == '__main__':
    m = Ch2()

    from mbuild.plot import Plot
    Plot(m, verbose=True, atoms=True, bonds=True, angles=False, dihedrals=False).show()


