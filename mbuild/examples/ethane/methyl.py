import numpy as np

from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import rotate_around_z, translate
from mbuild.testing.tools import get_fn


class Methyl(Compound):
    """A methyl group with one Port. """
    def __init__(self):
        super(Methyl, self).__init__(self)

        self.append_from_file(get_fn('methyl.mol2'))

        translate(self, -self.C[0])

        self.add(Port(anchor=self.C[0]), 'up')
        rotate_around_z(self.up, np.pi)
        translate(self.up, [0, -0.07, 0])

        self.add(Port(anchor=self.C[0]), 'down')
        translate(self.down, [0, -0.07, 0])

if __name__ == '__main__':
    methyl = Methyl()
    methyl.visualize(show_ports=True)

