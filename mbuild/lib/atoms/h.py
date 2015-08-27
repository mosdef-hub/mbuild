import numpy as np

import mbuild as mb


class H(mb.Compound):
    """A hydrogen atom with two overlayed ports."""
    def __init__(self):
        super(H, self).__init__()
        self.add(mb.Atom('H'), 'H')

        self.add(mb.Port(anchor=self.H), 'up')
        mb.rotate_around_z(self.up, np.pi)
        mb.translate(self.up, np.array([0, 0.07, 0]))

if __name__ == '__main__':
    m = H()
    m.visualize(show_ports=True)
