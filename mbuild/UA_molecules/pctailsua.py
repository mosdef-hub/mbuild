import mbuild as mb
import numpy as np

from mbuild.UA_molecules.ffaua import FFAUA
from mbuild.UA_molecules.alkyl_monomerua import AlkylMonomerUA

class PCTailsUA(mb.Compound):
    def __init__(self, tail_1_length, tail_2_length):
        super(PCTailsUA, self).__init__()
        
        self.add(mb.Polymer(AlkylMonomerUA(), 2), label='base')
        self.translate(-self['base'][1].pos)
        self['base']['up'].rotate(80*np.pi/180,[0,0,1])
        self['base']['up'].rotate(30*np.pi/180, [0,1,0])
        self['base'].add(mb.Port(anchor=self['base'][1],
            orientation=[(.143/2)*np.sin(19.5*np.pi/180),0,-(.143/2)*
                np.cos(19.5*np.pi/180)],
            separation=.143/2), label='side')
        
        self.add(FFAUA(tail_1_length), label='FFA[$]')
        mb.force_overlap(move_this=self['base'], 
                from_positions=self['base']['down'],
                to_positions=self['FFA'][0]['head']['down'])
        self.add(FFAUA(tail_2_length), label='FFA[$]')
        self['FFA'][1].translate(-self['FFA'][1]['head']['O'][1].pos)
        self['FFA'][1]['head']['down'].spin(np.pi,
                self['FFA'][1]['head']['down'].pos)
        mb.force_overlap(move_this=self['FFA'][1], 
                from_positions=self['FFA'][1]['head']['down'],
                to_positions=self['base']['side'])

if __name__ == '__main__':
    pctails = PCTailsUA(16,16)
    pctails.save('pctailsua.mol2', overwrite=True)
