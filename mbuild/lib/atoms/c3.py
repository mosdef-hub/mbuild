from __future__ import division

import numpy as np

import mbuild as mb


class C3(mb.Compound):
    """A tri-valent, planar carbon."""
    def __init__(self):
        super(C3, self).__init__()
        self.add(mb.Particle(name='C'))

        self.add(mb.Port(anchor=self[0]), 'up')
        self['up'].translate(np.array([0, 0.07, 0]))

        self.add(mb.Port(anchor=self[0]), 'down')
        self['down'].translate(np.array([0, 0.07, 0]))
        self['down'].spin(np.pi * 2/3, [0, 0, 1])

        self.add(mb.Port(anchor=self[0]), 'left')
        self['left'].translate(np.array([0, 0.07, 0]))
        self['left'].spin(-np.pi * 2/3, [0, 0, 1])


if __name__ == '__main__':
    m = C3()
    m.save('c3.mol2',overwrite=True)
