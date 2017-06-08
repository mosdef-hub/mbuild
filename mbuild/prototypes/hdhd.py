import mbuild as mb
import numpy as np

from mbuild.prototypes.ffa import FFA

class HDHD(mb.Compound):
    def __init__(self):
        super(HDHD,self).__init__()
        
        self.add(FFA(17), label='FFA[$]')
        self.add(FFA(17, hcap=True), label='FFA[$]')
        self['FFA'][1].remove(self['FFA'][1]['head']['OH'])
        self['FFA'][1].translate([0,.5,0])
        self['FFA'][0]['head']['down'].spin(np.pi, 
                self['FFA'][0]['head']['O'][1].pos-
                self['FFA'][0]['head']['C'].pos)
        self['FFA'][0].translate(-self['FFA'][0]['head']['O'][1].pos)
        self['FFA'][0]['head']['down'].rotate(25*np.pi/180, [0,0,1])
        
        mb.force_overlap(move_this=self['FFA'][1], 
                         from_positions=self['FFA'][1]['head']['down'], 
                         to_positions=self['FFA'][0]['head']['down'])
        
if __name__ == '__main__':
    hdhd = HDHD()
    hdhd.save('hdhd.mol2')

