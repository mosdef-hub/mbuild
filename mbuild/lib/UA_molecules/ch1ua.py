import mbuild as mb


class CH1UA(mb.Compound):
    """Creates a 3-port carbon for a united-atom model"""
    def __init__(self):
        super(CH1UA, self).__init__()
        self.add(mb.Particle(name='CH1'), label='C')
        self.add(mb.Port(anchor=self[0], orientation=[1, 1, 1], separation=.154/2), label='up')
        self.add(mb.Port(anchor=self[0], orientation=[-1, 1, -1], separation=.154/2), label='down')
        self.add(mb.Port(anchor=self[0], orientation=[-1, -1, 1], separation=.154/2), label='side')
