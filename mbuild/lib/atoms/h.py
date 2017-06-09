import numpy as np

import mbuild as mb


class H(mb.Compound):
    """A hydrogen atom with two overlayed ports."""
    def __init__(self):
        super(H, self).__init__()
        self.add(mb.Particle(name='H'))

        self.add(mb.Port(anchor=self[0]), 'up')
        self['up'].spin(np.pi, [0, 0, 1])
        self['up'].translate(np.array([0, 0.07, 0]))


if __name__ == '__main__':
    m = H()
    m.save('h.mol2', overwrite=True)
