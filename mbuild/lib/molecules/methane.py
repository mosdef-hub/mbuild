"""A methane molecule."""

import mbuild as mb


class Methane(mb.Compound):
    """A methane molecule."""

    def __init__(self):
        super(Methane, self).__init__()
        carbon = mb.Particle(name="C", element="C")
        self.add(carbon)

        hydrogen = mb.Particle(name="H", pos=[0.1, 0, -0.07], element="H")
        self.add(hydrogen)

        self.add_bond((self[0], self[1]), bond_order=1.0)

        self.add(mb.Particle(name="H", pos=[-0.1, 0, -0.07], element="H"))
        self.add(mb.Particle(name="H", pos=[0, 0.1, 0.07], element="H"))
        self.add(mb.Particle(name="H", pos=[0, -0.1, 0.07], element="H"))

        self.add_bond((self[0], self[2]), bond_order=1.0)
        self.add_bond((self[0], self[3]), bond_order=1.0)
        self.add_bond((self[0], self[4]), bond_order=1.0)
