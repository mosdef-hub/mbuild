import mbuild as mb
import numpy as np

from mbuild.UA_molecules.ch3ua import CH3UA
from mbuild.UA_molecules.ffaua import FFAUA
from mbuild.UA_molecules.alkyl_monomerua import AlkylMonomerUA
from mbuild.UA_molecules.alcua import ALCUA

class ISISUA(mb.Compound):
    def __init__(self):
        super(ISISUA, self).__init__()
        
        self.add(FFAUA(17), label='acid')
        self['acid'].name = 'acid'
        self['acid'].translate(-self['acid']['head']['O'][1].pos)
        self['acid']['head']['down'].rotate(125*np.pi/180, [1,0,0])
        self['acid']['head']['down'].spin(40*np.pi/180,
                self['acid']['head']['down'].pos)

        self.add(ALCUA(17), label='alcohol')
        self['alcohol'].name = 'alcohol'
        self['alcohol'].translate(-self['alcohol']['tail'][0].pos)
        self['alcohol'].add(mb.Port(anchor=self['alcohol']['tail'][0], 
            orientation=self['alcohol']['head'][0].pos,
            separation=.15/2), label='up')
        self['alcohol']['up'].spin(-3*np.pi/4, self['alcohol']['up'].pos)
        self['alcohol'].remove(self['alcohol']['head'])
        mb.force_overlap(move_this=self['alcohol'], 
                from_positions=self['alcohol']['up'],
                to_positions=self['acid']['head']['down'])
        
        self['alcohol'].add(mb.Port(anchor=self['alcohol'][16], 
            orientation=[-1,0.8,-1], 
            separation=.14/2), label='down')
        self['alcohol'].add(CH3UA(), label='methyl')
        mb.force_overlap(move_this=self['alcohol']['methyl'], 
                from_positions=self['alcohol']['methyl']['up'],
                to_positions=self['alcohol']['down'])
        
        self.translate(-self['acid']['tail'][0].pos)
        self['acid']['tailcap'].rotate(np.pi, self['acid']['head']['C'].pos)
        self['acid']['tail'].rotate(np.pi, self['acid']['head']['C'].pos)
        
        self.add(mb.Port(anchor=self['acid'][17], 
            orientation=[3,-0.5,-3],
            separation=.14/2), label='acid_branch_port')
        self.add(CH3UA(), label='acid_branch')
        mb.force_overlap(move_this=self['acid_branch'], 
                from_positions=self['acid_branch']['up'],
                to_positions=self['acid_branch_port'])
        
        mb.z_axis_transform(self, new_origin = self['acid']['head']['O'][0], 
                            point_on_z_axis = self['acid']['tailcap'][0])
        self.rotate(np.pi/2, [0,0,1])
        self.rotate(np.pi, [0,1,0])
        self.name = 'ISIS'

if __name__ == '__main__':
    isis = ISISUA()
    isis.save('isisua.mol2', overwrite=True)
