import mbuild as mb
import numpy as np

class OH(mb.Compound):
    """Creates a hydroxyl group"""
    def __init__(self):
        super(OH,self).__init__()
        self.add(mb.Particle(name='O', pos = [0,0,0]), label='O')
        self.add(mb.Particle(name='H', pos = [.096, 0, 0]), label='H')
        self.add_bond((self[0], self[1]))

        self.add(mb.Port(anchor=self[0]), label='cbond')
        self['cbond'].translate([-.143/2, 0, 0])
        theta = (180 - 104.5) * np.pi / 180
        self['cbond'].rotate(theta, around=[0,1,0])

if __name__ == '__main__':
    oh = OH()
    oh.energy_minimization()
    oh.visualize()
    oh.save('oh.mol2')
