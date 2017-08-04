import mbuild as mb
import numpy as np

from mbuild.lib.prototypes.pctails import PCTails
from mbuild.lib.prototypes.pchead import PCHead


class DSPC(mb.Compound):
    def __init__(self):
        super(DSPC, self).__init__()
        
        self.add(PCHead(), label='headgroup')
        self.add(PCTails(18, 18), label='ffatails')
        
        mb.force_overlap(move_this=self['ffatails'], from_positions=self['ffatails']['base']['up'],
                         to_positions=self['headgroup']['alkyl_split']['up'])
        mb.z_axis_transform(self, self['headgroup']['N4']['N'][0], self['headgroup']['PO4']['P'][0])
        self.spin(np.pi, [1, 0, 0])

if __name__ == '__main__':
    dspc = DSPC()
    dspc.save('dspc.mol2', overwrite=True)
