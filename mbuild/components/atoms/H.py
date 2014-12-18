import numpy as np

from mbuild.atom import Atom
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import translate, rotate_around_z


class H(Compound):
    """A hydrogen atom with two overlayed ports. """
    def __init__(self):
        super(H, self).__init__(self)
        self.add(Atom('H'), 'H')

        self.add(Port(anchor=self.H), 'up')
        rotate_around_z(self.up, np.pi)
        translate(self.up, np.array([0, 0.07, 0]))

        self.add(Port(anchor=self.H), 'down')
        translate(self.down, np.array([0, 0.07, 0]))


if __name__ == '__main__':
    m = H()
    m.visualize(show_ports=True)
