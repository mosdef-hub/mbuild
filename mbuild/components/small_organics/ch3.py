import numpy as np

from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import translate, rotate_around_z


class Ch3(Compound):
    """A methyl group. """
    def __init__(self):
        super(Ch3, self).__init__(self)

        self.append_from_file('ch3.pdb')
        translate(self, -self.C[0])

        self.add(Port(anchor=self.C[0]), 'up')
        rotate_around_z(self.up, np.pi)
        translate(self.up, [0, 0.07, 0])

        self.add(Port(anchor=self.C[0]), 'down')
        translate(self.down, [0, 0.07, 0])


if __name__ == '__main__':
    m = Ch3()
    m.visualize(show_ports=True)

