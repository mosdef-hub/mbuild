import mbuild as mb
import numpy as np


class OH(mb.Compound):
    """Creates a hydroxyl group"""
    def __init__(self):
        super(OH, self).__init__()
        self.add(mb.Particle(name='OA', pos=[0, 0, 0]), label='O')
        self.add(mb.Particle(name='H', pos=[.096, 0, 0]), label='H')
        self.add_bond((self[0], self[1]))

        self.add(mb.Port(anchor=self[0], orientation=[-1, np.tan(72*np.pi/180), 0], separation=.135/2), label='up')

if __name__ == '__main__':
    oh = OH()
    oh.energy_minimization()
    oh.save('oh.mol2', overwrite=True)
