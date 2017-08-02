import mbuild as mb
import numpy as np

from mbuild.lib.UA_molecules.ffaua import FFAUA
from mbuild.lib.UA_molecules.alkyl_monomerua import AlkylMonomerUA
from mbuild.lib.UA_molecules.ch1ua import CH1UA


class PCTailsUA(mb.Compound):
    def __init__(self, tail_1_length, tail_2_length):
        super(PCTailsUA, self).__init__()
        
        self.add(CH1UA(), label='CH1')
        self.add(FFAUA(tail_1_length), label='FFA[$]')
        self['FFA'][0].translate(-self['FFA'][0]['head']['O'][0].pos)
        self['FFA'][0]['head']['down'].spin(-np.pi/2, self['FFA'][0]['head']['down'].pos)
        mb.force_overlap(move_this=self['CH1'], from_positions=self['CH1']['up'],
                         to_positions=self['FFA'][0]['head']['down'])
        self.add(AlkylMonomerUA(), label='CH2')
        self['CH2']['up'].spin(-np.pi/2, self['CH2']['up'].pos)
        self['CH2']['down'].spin(-np.pi/2, self['CH2']['down'].pos)
        mb.force_overlap(move_this=self['CH2'], from_positions=self['CH2']['up'],
                         to_positions=self['CH1']['down'])
    
        self.add(FFAUA(tail_2_length), label='FFA[$]')
        self['FFA'][1].translate(-self['FFA'][1]['head']['O'][1].pos)
        self['FFA'][1]['head']['down'].spin(np.pi, self['FFA'][1]['head']['down'].pos)
        mb.force_overlap(move_this=self['FFA'][1], from_positions=self['FFA'][1]['head']['down'],
                         to_positions=self['CH2']['down'])
        self['FFA'][0].name = 'ffatail'
        self['FFA'][1].name = 'ffatail'

if __name__ == '__main__':
    pctailsua = PCTailsUA(14, 14)
    pctailsua.save('pctailsua.mol2', overwrite=True)
