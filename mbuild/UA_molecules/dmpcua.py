import mbuild as mb
import numpy as np

from mbuild.UA_molecules.pctailsua import PCTailsUA
from mbuild.UA_molecules.pcheadua import PCHeadUA

class DMPCUA(mb.Compound):
    def __init__(self):
        super(DMPCUA, self).__init__()
        
        self.add(PCHeadUA(), label='headgroup')
        self.add(PCTailsUA(14,14), label='ffatails')
        
        mb.force_overlap(move_this=self['ffatails'],
                        from_positions=self['ffatails']['base']['up'],
                        to_positions=self['headgroup']['alkyl_split']['up'])
    
        mb.z_axis_transform(self, new_origin=self['headgroup'][0],
                point_on_z_axis=self['ffatails']['FFA'][0][0])
        self.rotate(np.pi, [1,0,0])

if __name__ == '__main__':
    dmpc = DMPCUA()
    dmpc.save('dmpcua.mol2', overwrite=True)

