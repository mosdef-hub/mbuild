import mbuild as mb
import numpy as np

from mbuild.UA_molecules.pctailsua  import PCTailsUA
from mbuild.UA_molecules.pcheadua import PCHeadUA

class DSPCUA(mb.Compound):
    def __init__(self):
        super(DSPCUA, self).__init__()
        
        self.add(PCHeadUA(), label='headgroup')
        self.add(PCTailsUA(18,18), label='ffatails')
        
        mb.force_overlap(move_this=self['ffatails'],
                        from_positions=self['ffatails']['CH1']['side'],
                        to_positions=self['headgroup']['alkyl_split']['up'])

        mb.z_axis_transform(self, new_origin=self['headgroup']['N'],
                point_on_z_axis=self['ffatails']['FFA'][1][12],
                point_on_zx_plane=self['ffatails']['FFA'][0]['C'])
        self.rotate(np.pi, [1,0,0])
        self.name = 'DSPC'

if __name__ == '__main__':
    dspc = DSPCUA()
    dspc.save('dspcua.mol2', overwrite=True)

