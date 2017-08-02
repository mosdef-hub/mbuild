import mbuild as mb
import numpy as np


class Phosphate(mb.Compound):
    """Creates a phosphate ion with 2 ports attached"""
    def __init__(self):
        super(Phosphate, self).__init__()
        
        self.add(mb.Particle(name='P'), label='P')
        self.add(mb.Particle(name='O', pos=[0, .16, 0]), label='O[$]')
        self.add(mb.Particle(name='O', pos=[0, -.16, 0]), label='O[$]')
        self.add(mb.Particle(name='O', pos=[.16, 0, 0]), label='O[$]')
        self.add(mb.Particle(name='O', pos=[-.16, 0, 0]), label='O[$]')
        self.add_bond((self['P'], self['O[0]']))
        self.add_bond((self['P'], self['O[1]']))
        self.add_bond((self['P'], self['O[2]']))
        self.add_bond((self['P'], self['O[3]']))
        
        alpha = (180-109.5)*np.pi/180
        
        self['O'][1].rotate(alpha, [1, 0, 0])
        self['O'][2].rotate(-np.pi/6, [0, 1, 0])
        self['O'][3].rotate(np.pi/6, [0, 1, 0])
        
        self['O'][2].rotate(19.5*np.pi/180, [self['O'][2].pos[0], 0, -self['O'][2].pos[2]])
        self['O'][3].rotate(-19.5*np.pi/180, [self['O'][3].pos[0], 0, -self['O'][3].pos[2]])
        
        mb.x_axis_transform(self, new_origin=self['O'][1], point_on_x_axis=self['O'][0],
                            point_on_xy_plane=self['P'])
        
        self.add(mb.Port(anchor=self['O'][0], orientation=[1, np.tan(36.75*np.pi/180), 0],
                         separation=.143/2), label='up')
        self.add(mb.Port(anchor=self['O'][1], orientation=[-1, np.tan(36.75*np.pi/180), 0],
                         separation=.143/2), label='down')

        # follow .itp naming convention
        self['O'][0].name = 'OA'
        self['O'][1].name = 'OA'
        self['O'][2].name = 'OM'
        self['O'][3].name = 'OM'

if __name__ == '__main__':
    po4 = Phosphate()
    po4.save('po4.mol2', overwrite=True)
