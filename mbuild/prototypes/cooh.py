import mbuild as mb
import numpy as np

from mbuild.prototypes.OH import OH

class COOH(mb.Compound):
    """Creates headgroup of a carboxylic acid"""
    def __init__(self):
        super(COOH,self).__init__()
        self.add(mb.Particle(name='C'), label='C')
        self.add(mb.Particle(name='O', pos = [.123, 0, 0]),label='CO')
        self.add_bond((self[0],self[1]))
                                                        
        self.add(mb.Port(anchor=self[0]), label='up')
        self['up'].translate([-.132/2,0,0])
        theta = (111 / 2) * np.pi / 180
        self['up'].rotate(theta, around=[0,0,1])
        
        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([-.15/2,0,0])
        self['down'].rotate(-theta, around=[0,0,1])
        self.add(OH(),label='hydroxyl')
        mb.force_overlap(move_this=self['hydroxyl'],
                from_positions=self['hydroxyl']['cbond'],
                to_positions=self['up'])

if __name__ == '__main__':
    cooh = COOH()
    cooh.energy_minimization()
    cooh.save('cooh.mol2')
