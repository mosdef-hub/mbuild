import mbuild as mb
import numpy as np

class CH1UA(mb.Compound):
    """Creates a 3-port carbon for a united-atom model"""
    def __init__(self):
        super(CH1UA, self).__init__()
        self.add(mb.Particle(name='CH1', pos=[0,0,0]), label='C[$]')
        
        theta = (180-109.5)/2 * np.pi/180
        self.add(mb.Port(anchor=self[0]), label='up')
        self['up'].translate([0, -0.154/2, 0])
        self['up'].rotate(theta, around=[1,0,0])
        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, 0.154/2, 0])
        self['down'].rotate(np.pi, around=[0,1,0])
        self['down'].rotate(-theta, around=[1,0,0])
        self.add(mb.Port(anchor=self[0], orientation=[-1,0,0],
            separation = .154/2), label='side')
        self['side'].rotate(theta, around=[0,1,0])
        self['side'].spin(np.pi, around=self['side'].pos)

