__author__ = 'sallai'
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import *
from mbuild.testing.tools import get_fn


class Methyl(Compound):
    """ """
    def __init__(self):
        Compound.__init__(self)

        self.append_from_file(get_fn('methyl.pdb'))

        translate(self, -self.C[0])

        self.add(Port(anchor=self.C[0]), 'up')
        rotate_around_z(self.up, np.pi)
        translate(self.up, np.array([0, -0.07, 0]))

        self.add(Port(anchor=self.C[0]), 'down')
        translate(self.down, np.array([0, -0.07, 0]))

if __name__ == '__main__':
    methyl = Methyl()

    methyl.visualize()

