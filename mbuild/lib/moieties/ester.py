import numpy as np

import mbuild as mb


class Ester(mb.Compound):
    """A ester group -C(=O)O-. """
    def __init__(self):
        super(Ester, self).__init__()

        mb.load('ester.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[2]), 'up')
        self['up'].spin(np.pi / 2, [0, 0, 1])
        self['up'].translate_to(np.array([0.07, 0, 0]))

        self.add(mb.Port(anchor=self[0]), 'down')
        self['down'].spin(np.pi / 2, [0, 0, 1])
        self['down'].translate(np.array([-0.07, 0, 0]))

if __name__ == '__main__':
    m = Ester()
    m.save('ester.mol2', overwrite=True)
