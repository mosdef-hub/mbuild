import mbuild as mb
import numpy as np

from mbuild.lib.UA_molecules.pctailsua import PCTailsUA
from mbuild.lib.UA_molecules.pcheadua import PCHeadUA


class DSPCUA(mb.Compound):
    def __init__(self):
        super(DSPCUA, self).__init__()
        
        self.add(PCHeadUA(), label='headgroup')
        self.add(PCTailsUA(18, 18), label='ffatails')
        
        mb.force_overlap(move_this=self['ffatails'],
                         from_positions=self['ffatails']['CH1']['side'],
                         to_positions=self['headgroup']['alkyl_split']['up'])
        self.translate(-self['ffatails']['CH1']['C'].pos)
        self['headgroup'].rotate(-60*np.pi/180, [0, 1, 0])

        mb.z_axis_transform(self, new_origin=self['ffatails']['FFA'][1][5],
                            point_on_z_axis=self['ffatails']['FFA'][1][7],
                            point_on_zx_plane=self['ffatails']['FFA'][1][6])
        self.rotate(np.pi, [1, 0, 0])
        self.name = 'DSPC'

if __name__ == '__main__':
    dspc = DSPCUA()
    dspc.save('dspcua.mol2', overwrite=True)
