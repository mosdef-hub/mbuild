import mbuild as mb
import numpy as np


class AlkylMonomer(mb.Compound):
    """Creates a monomer for an alkyl chain"""
    def __init__(self):
        super(AlkylMonomer, self).__init__()
        self.add(mb.Particle(name='C', pos=[0, 0, 0]), label='C[$]')
        self.add(mb.Particle(name='H', pos=[-0.109, 0, 0]), label='HC[$]')
        self.add(mb.Particle(name='H', pos=[0.109, 0, 0]), label='HC[$]')
        theta = 0.5 * (180 - 109.5) * np.pi / 180
        self['HC'][0].rotate(theta, around=[0, 1, 0])
        self['HC'][1].rotate(-theta, around=[0, 1, 0])
        self.add_bond((self[0], self['HC'][0]))
        self.add_bond((self[0], self['HC'][1]))

        self.add(mb.Port(anchor=self[0]), label='up')
        self['up'].translate([0, -0.154/2, 0])
        self['up'].rotate(theta, around=[1, 0, 0])
        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, 0.154/2, 0])
        self['down'].rotate(np.pi, around=[0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])

if __name__ == '__main__':
    ch2 = AlkylMonomer()
    ch2.save('alkylmonomer.mol2')
