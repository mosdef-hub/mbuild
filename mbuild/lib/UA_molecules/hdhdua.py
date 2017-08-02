import mbuild as mb
import numpy as np

from mbuild.lib.UA_molecules import CH3UA, AlkylMonomerUA
from mbuild.lib.prototypes import COOH


class HDHDUA(mb.Compound):
    def __init__(self):
        super(HDHDUA, self).__init__()
        
        self.add(CH3UA(), label='tail_end[$]')
        self.add(mb.Polymer(AlkylMonomerUA(), 16), label='tail[$]')
        self.add(COOH(ester=True, cis=True), label='head')
        self['head']['up'].spin(-np.pi/2, self['head']['up'].pos)
        self.add(mb.Polymer(AlkylMonomerUA(), 15), label='tail[$]')
        self.add(CH3UA(), label='tail_end[$]')
        
        mb.force_overlap(move_this=self['tail'][0], from_positions=self['tail'][0]['up'],
                         to_positions=self['head']['down'])
        mb.force_overlap(move_this=self['tail'][1], from_positions=self['tail'][1]['down'],
                         to_positions=self['head']['up'])
      
        mb.force_overlap(move_this=self['tail_end'][0], from_positions=self['tail_end'][0]['up'],
                         to_positions=self['tail'][0]['up'])
        mb.force_overlap(move_this=self['tail_end'][1], from_positions=self['tail_end'][1]['up'],
                         to_positions=self['tail'][1]['up'])
        
        mb.z_axis_transform(compound=self, new_origin=self['head']['C'], point_on_z_axis=self['tail'][1][3],
                            point_on_zx_plane=self['tail'][1][1])
        self.spin(np.pi, [1, 0, 0])
        self.name = 'HDHD'

if __name__ == '__main__':
    hdhd = HDHDUA()
    hdhd.save('hdhdua.mol2', overwrite=True)
