import numpy as np

import mbuild as mb


class Ester(mb.Compound):
    """A ester group -C(=O)O-. """
    def __init__(self):
        super(Ester, self).__init__()

        mb.load('ester.pdb', compound=self, relative_to_module=self.__module__)
        mb.translate(self, -self[0].pos)

        self.add(mb.Port(anchor=self[2]), 'up')
        mb.rotate_around_z(self['up'], np.pi / 2)
        mb.translate_to(self['up'], self[2].pos + np.array([0.07, 0, 0]))

        self.add(mb.Port(anchor=self[0]), 'down')
        mb.rotate_around_z(self['down'], np.pi / 2)
        mb.translate(self['down'], np.array([-0.07, 0, 0]))

if __name__ == '__main__':
    m = Ester()
    m.visualize(show_ports=True)

