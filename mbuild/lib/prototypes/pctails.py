import mbuild as mb
import numpy as np

from mbuild.lib.prototypes.ffa import FFA
from mbuild.lib.prototypes.alkyl_monomer import AlkylMonomer


class PCTails(mb.Compound):
    """Creates acyl tailgroups for a phosphatidylcholine molecule"""
    def __init__(self, tail_1_length, tail_2_length):
        super(PCTails, self).__init__()
        
        self.add(mb.Polymer(AlkylMonomer(), 2), label='base')
        self.translate(-self['base'][3].pos)
        self['base'].add(mb.Port(anchor=self['base'][3], orientation=self['base'][5].pos,
                                 separation=.143/2), label='side')
        self.remove(self['base'][5])
        self['base']['side'].spin(np.pi, self['base']['side'].pos)
        self.translate(-self['base'][0].pos)
        self['base']['down'].spin(np.pi/2, self['base']['down'].pos)
        
        self.add(FFA(tail_1_length, ester=True), label='FFA[$]')
        mb.force_overlap(move_this=self['base'], from_positions=self['base']['down'],
                         to_positions=self['FFA'][0]['head']['down'])
        self.add(FFA(tail_2_length, ester=True), label='FFA[$]')
        mb.force_overlap(move_this=self['FFA'][1], from_positions=self['FFA'][1]['head']['down'],
                         to_positions=self['base']['side'])

if __name__ == '__main__':
    pctails = PCTails(16, 16)
    pctails.save('pctails.mol2', overwrite=True)
