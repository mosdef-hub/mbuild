import mbuild as mb
import numpy as np

from mbuild.UA_molecules.pctailsua import PCTailsUA
from mbuild.UA_molecules.pcheadua import PCHeadUA

class DPPCUA(mb.Compound):
    def __init__(self):
        super(DPPCUA, self).__init__()
        
        self.add(PCHeadUA(), label='headgroup')
        self.add(PCTailsUA(16,16), label='ffatails')
        
        mb.force_overlap(move_this=self['ffatails'],
                        from_positions=self['ffatails']['CH1']['side'],
                        to_positions=self['headgroup']['alkyl_split']['up'])

        mb.z_axis_transform(self, new_origin=self['ffatails']['CH1'],
                point_on_z_axis=self['ffatails']['FFA'][0][4],
                point_on_zx_plane=self['ffatails']['FFA'][1]['C'])
        self.rotate(np.pi, [1,0,0])
        self.name = 'DPPC'
    
if __name__ == '__main__':
    dppc = DPPCUA()
    dppc.save('dppcua.mol2', overwrite=True)

