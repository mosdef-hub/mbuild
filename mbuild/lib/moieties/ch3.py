import mbuild as mb


class CH3(mb.Compound):
    """A methyl group. """
    def __init__(self):
        super(CH3, self).__init__()

        mb.load('ch3.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)  # Move carbon to origin.

        self.add(mb.Port(anchor=self[0]), 'up')
        self['up'].translate([0, -0.07, 0])

if __name__ == '__main__':
    m = CH3()
    m.save('ch3.mol2', overwrite=True)
