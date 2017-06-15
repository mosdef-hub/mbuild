__author__ = 'jonestj1'

import mbuild as mb

class PegMonomer(mb.Compound):
    def __init__(self):
        super(PegMonomer, self).__init__()

        mb.load('peg_monomer.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0]), 'down')
        self['down'].translate([0, -0.07, 0])

        self.add(mb.Port(anchor=self[6]), 'up')
        self['up'].translate([0, 0.073, 0])

if __name__ == '__main__':
    peg = PegMonomer()
    peg.save('peg.mol2')
