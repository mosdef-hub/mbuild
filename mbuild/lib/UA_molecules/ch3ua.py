import mbuild as mb


class CH3UA(mb.Compound):
    """A methyl group for united-atom model. """
    def __init__(self):
        super(CH3UA, self).__init__()
        self.add(mb.Particle(name='C'), label='C')
        self.add(mb.Port(anchor=self[0]), 'up')
        self['up'].translate([0, -0.07, 0])

        self[0].name = 'CH3'

if __name__ == '__main__':
    m = CH3UA()
