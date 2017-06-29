import mbuild as mb

class CH3UA(mb.Compound):
    """A methyl group for united-atom model. """
    def __init__(self):
        super(CH3UA, self).__init__()

        #mb.load('ch3ua.pdb', compound=self, relative_to_module=self.__module__)
        #self.translate(-self[0].pos)  # Move carbon to origin.

        self.add(mb.Particle(name='C'), label='C')
        self.add(mb.Port(anchor=self[0]), 'up')
        self['up'].translate([0, -0.07, 0])

        #self.remove([self[1], self[2], self[3]])

        self[0].name = 'CH3'
        

if __name__ == '__main__':
    m = CH3UA()

