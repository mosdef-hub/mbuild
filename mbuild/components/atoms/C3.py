from __future__ import division

import numpy as np

from mbuild.atom import Atom
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import translate, rotate_around_z, translate_to, \
    _extract_atom_positions


class C3(Compound):
    """A tri-valent, planar carbon. """
    def __init__(self):
        super(C3, self).__init__(self)
        self.add(Atom('C'), 'C')

        self.add(Port(anchor=self.C), 'up')
        translate(self.up, np.array([0, 0.07, 0]))

        self.add(Port(anchor=self.C), 'down')
        translate(self.down, np.array([0, 0.07, 0]))
        rotate_around_z(self.down, np.pi * 2/3)

        self.add(Port(anchor=self.C), 'left')
        translate(self.left, np.array([0, 0.07, 0]))
        rotate_around_z(self.left, -np.pi * 2/3)


if __name__ == '__main__':
    m = C3()
    m.visualize(show_ports=True)
