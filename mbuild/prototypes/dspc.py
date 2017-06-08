import mbuild as mb
import numpy as np

from mbuild.prototypes.pctails import PCTails
from mbuild.prototypes.pchead import PCHead

class DSPC(mb.Compound):
    def __init__(self):
        super(DSPC, self).__init__()
        
        self.add(PCHead(), label='headgroup')
        self.add(PCTails(18,18), label='ffatails')
        
        mb.force_overlap(move_this=self['ffatails'],
                        from_positions=self['ffatails']['base']['up'],
                        to_positions=self['headgroup']['alkyl_split']['up'])
    
if __name__ == '__main__':
    dspc = DSPC()
    dspc.save('dspc.mol2', overwrite=True)

