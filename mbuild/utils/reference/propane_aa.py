import mbuild as mb

class Propane_aa(mb.Compound):
    """A propane molecule with one hydrogen omitted for bonding """
    def __init__(self):
        super(Propane_aa, self).__init__()

        mb.load('propane_aa.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)  # Move carbon to origin.

        self.add(mb.Port(anchor=self[0]), 'up')
        self['up'].translate([0, 0.07, 0])

if __name__ == '__main__':
    m = Propane_aa()
    m.save('propane_aa.mol2', overwrite=True)
