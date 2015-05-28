from __future__ import division

import numpy as np

import mbuild as mb


class C3(mb.Compound):
    """A tri-valent, planar carbon."""
    def __init__(self):
        super(C3, self).__init__()
        self.add(mb.Atom('C'), 'C')

        self.add(mb.Port(anchor=self.C), 'up')
        mb.translate(self.up, np.array([0, 0.07, 0]))

        self.add(mb.Port(anchor=self.C), 'down')
        mb.translate(self.down, np.array([0, 0.07, 0]))
        mb.rotate_around_z(self.down, np.pi * 2/3)

        self.add(mb.Port(anchor=self.C), 'left')
        mb.translate(self.left, np.array([0, 0.07, 0]))
        mb.rotate_around_z(self.left, -np.pi * 2/3)


if __name__ == '__main__':
    m = C3()
    m.visualize(show_ports=True)
