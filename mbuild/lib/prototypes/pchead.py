import mbuild as mb
import numpy as np

from mbuild.examples.methane.methane import Methane
from mbuild.lib.atoms.n4 import N4
from mbuild.lib.prototypes.po4 import Phosphate
from mbuild.lib.prototypes.alkyl_monomer import AlkylMonomer


class PCHead(mb.Compound):
    def __init__(self):
        super(PCHead, self).__init__()
        
        self.add(N4(), label='N4')
        
        for i in range(0, 3):
            self.add(Methane(), label='methyl[$]')
            self['methyl'][i].add(mb.Port(anchor=self['methyl'][i][0], orientation=self['methyl'][i][1].pos,
                                          separation=.077), label='up')
            self['methyl'][i].remove(self['methyl'][i][1])
            self['methyl'][i].translate([-.5, 0, 0])
        
            mb.force_overlap(move_this=self['methyl'][i], from_positions=self['methyl'][i]['up'],
                             to_positions=self['N4']['port_'+str(i)])
        
        self.add(mb.Polymer(AlkylMonomer(), 2), label='alkyl_body')
        mb.force_overlap(move_this=self['alkyl_body'], 
                         from_positions=self['alkyl_body']['down'],
                         to_positions=self['N4']['port_3'])
        self.energy_minimization()
    
        self.add(Phosphate(), label='PO4')
        mb.force_overlap(move_this=self['PO4'], from_positions=self['PO4']['down'],
                         to_positions=self['alkyl_body']['up'])
        
        self.add(AlkylMonomer(), label='alkyl_split')
        self['alkyl_split']['down'].spin(np.pi/2, self['alkyl_split']['down'].pos)
        mb.force_overlap(move_this=self['alkyl_split'], from_positions=self['alkyl_split']['down'],
                         to_positions=self['PO4']['up'])
        
if __name__ == '__main__':
    pchead = PCHead()
    pchead.save('pchead.mol2', overwrite=True)
