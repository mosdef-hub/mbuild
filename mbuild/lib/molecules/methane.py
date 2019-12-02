# -*- coding: utf-8 -*-


# -- ==methane== --
import mbuild as mb


class Methane(mb.Compound):
    def __init__(self):
        super(Methane, self).__init__()
        carbon = mb.Particle(name='C')
        self.add(carbon)

        hydrogen = mb.Particle(name='H', pos=[0.1, 0, -0.07])
        self.add(hydrogen)

        self.add_bond((self[0], self[1]))

        self.add(mb.Particle(name='H', pos=[-0.1, 0, -0.07]))
        self.add(mb.Particle(name='H', pos=[0, 0.1, 0.07]))
        self.add(mb.Particle(name='H', pos=[0, -0.1, 0.07]))

        self.add_bond((self[0], self[2]))
        self.add_bond((self[0], self[3]))
        self.add_bond((self[0], self[4]))

# -- ==methane== --