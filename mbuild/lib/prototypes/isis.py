import mbuild as mb
import numpy as np

from mbuild.lib.moieties.ch3 import CH3
from mbuild.lib.prototypes.ffa import FFA
from mbuild.lib.prototypes.alc import ALC


class ISIS(mb.Compound):
    def __init__(self):
        super(ISIS, self).__init__()
        
        self.add(FFA(17), label='acid')
        self['acid'].translate(-self['acid']['head']['O'][1].pos)
        self['acid']['head']['down'].rotate(125*np.pi/180, [1, 0, 0])
        self.add(ALC(17), label='alcohol')
        self['alcohol'].translate(-self['alcohol']['tail'][0].pos)
        self['alcohol'].add(mb.Port(anchor=self['alcohol']['tail'][0], orientation=self['alcohol']['head'][0].pos,
                                    separation=.15/2), label='up')
        self['alcohol']['up'].spin(-3*np.pi/4, self['alcohol']['up'].pos)
        self['alcohol'].remove(self['alcohol']['head'])
        mb.force_overlap(move_this=self['alcohol'], from_positions=self['alcohol']['up'],
                         to_positions=self['acid']['head']['down'])
        
        self['alcohol'].add(mb.Port(anchor=self['alcohol']['tail'][45],
                                    orientation=(self['alcohol']['tail'][46].pos - self['alcohol']['tail'][45].pos),
                                    separation=.14/2), label='down')
        self['alcohol'].remove(self['alcohol']['tail'][46])
        self['alcohol'].add(CH3(), label='methyl')
        mb.force_overlap(move_this=self['alcohol']['methyl'], from_positions=self['alcohol']['methyl']['up'],
                         to_positions=self['alcohol']['down'])
        
        self.translate(-self['acid']['tail'][0].pos)
        self['acid']['tailcap'].rotate(np.pi, self['acid']['head']['C'].pos)
        self['acid']['tail'].rotate(np.pi, self['acid']['head']['C'].pos)
        
        self.add(mb.Port(anchor=self['acid']['tail'][42], 
                         orientation=(self['acid']['tail'][44].pos - self['acid']['tail'][42].pos),
                         separation=.14/2), label='acid_branch_port')
        self['acid'].remove(self['acid']['tail'][44])
        self.add(CH3(), label='acid_branch')
        mb.force_overlap(move_this=self['acid_branch'], from_positions=self['acid_branch']['up'],
                         to_positions=self['acid_branch_port'])
        
        mb.z_axis_transform(self, new_origin=self['acid']['head']['O'][0],
                            point_on_z_axis=self['acid']['tailcap'][0])
        self.rotate(np.pi/2, [0, 0, 1])
        self.rotate(np.pi, [0, 1, 0])

if __name__ == '__main__':
    isis = ISIS()
    isis.save('isis.mol2', overwrite=True)
