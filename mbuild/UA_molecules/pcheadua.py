import mbuild as mb
import numpy as np

from mbuild.lib.atoms.n4 import N4
from mbuild.UA_molecules.po4ua import Phosphate
from mbuild.UA_molecules.alkyl_monomerua import AlkylMonomerUA
from mbuild.UA_molecules.ch3ua import CH3UA

class PCHeadUA(mb.Compound):
    def __init__(self):
        super(PCHeadUA,self).__init__()
        
        self.add(N4(), label='N4')
        
        for i in range(0,3):
            self.add(CH3UA(), label='methyl[$]')
        
            mb.force_overlap(move_this=self['methyl'][i],
                    from_positions=self['methyl'][i]['up'],
                    to_positions=self['N4']['port_'+str(i)])
        
        self.add(mb.Polymer(AlkylMonomerUA(),2), label='alkyl_body')
        mb.force_overlap(move_this=self['alkyl_body'], 
                         from_positions=self['alkyl_body']['down'],
                         to_positions=self['N4']['port_3'])
    
        self.add(Phosphate(), label='PO4')
        #self['PO4']['down'].spin(np.pi/2, self['PO4']['down'].pos)   
        #single P-O bond behavior
        mb.force_overlap(move_this=self['PO4'], 
                from_positions=self['PO4']['down'],
                to_positions=self['alkyl_body']['up'])
        
        self.add(AlkylMonomerUA(), label='alkyl_split')
        self['alkyl_split']['down'].spin(np.pi/2,
                self['alkyl_split']['down'].pos)
        mb.force_overlap(move_this=self['alkyl_split'], 
                         from_positions=self['alkyl_split']['down'],
                         to_positions=self['PO4']['up'])
        
if __name__ == '__main__':
    pchead=PCHeadUA()
    pchead.save('pcheadua.mol2', overwrite=True)
