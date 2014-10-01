__author__ = 'sallai'
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import *
from mbuild.testing.tools import get_fn

class Ch2(Compound):
    """ """
    def __init__(self):
        Compound.__init__(self)

        # Look for data file in same directory as this python module.
        self.append_from_file(get_fn('ch2.pdb'))

        self.add(Port(anchor=self.C[0]), 'up')
        translate(self.up, np.array([0, 0.07, 0]))

        self.add(Port(anchor=self.C[0]), 'down')
        translate(self.down, np.array([0, -0.07, 0]))

if __name__ == '__main__':
    m = Ch2()

    m.visualize()

