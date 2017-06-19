import mbuild as mb
import numpy as np

from mbuild.UA_molecules.ffaua import FFAUA
from mbuild.UA_molecules.alcua import ALCUA

class HDHDUA(mb.Compound):
    def __init__(self):
        super(HDHDUA, self).__init__()
        
        self.add(FFAUA(17), label='acid')
        self['acid'].translate(-self['acid']['head']['O'][1].pos)
        self['acid']['head']['down'].rotate(125*np.pi/180, [1,0,0])
        self['acid']['head']['down'].spin(40*np.pi/180, 
                self['acid']['head']['down'].pos)
        
        self.add(ALCUA(17), label='alcohol')
        self['alcohol'].translate(-self['alcohol']['tail'][0].pos)
        self['alcohol'].add(mb.Port(anchor=self['alcohol']['tail'][0], 
            orientation=self['alcohol']['head'][0].pos,
            separation=.15/2), label='up')
        self['alcohol']['up'].spin(-3*np.pi/4, self['alcohol']['up'].pos)
        self['alcohol'].remove(self['alcohol']['head'])
        mb.force_overlap(move_this=self['alcohol'], 
                from_positions=self['alcohol']['up'],
                to_positions=self['acid']['head']['down'])
        
        self.translate(-self['acid']['tail'][0].pos)
        self['acid']['tailcap'].rotate(np.pi, self['acid']['head']['C'].pos)
        
        self['acid']['tail'].rotate(np.pi, self['acid']['head']['C'].pos)
        
        mb.z_axis_transform(self, new_origin = self['acid']['head']['O'][0], 
                           point_on_z_axis = self['acid']['tailcap'][0])
        self.rotate(np.pi/2, [0,0,1])
        self.rotate(np.pi, [0,1,0])
        
if __name__ == '__main__':
    hdhd = HDHDUA()
    hdhd.save('hdhdua.mol2', overwrite=True)
