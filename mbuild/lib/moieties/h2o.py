import mbuild as mb


class H2O(mb.Compound):
    """A water molecule. """
    def __init__(self):
        super(H2O, self).__init__()

        self.add(mb.Particle(name='O', pos=[1.0203, 0.7604, 1.2673]), 'OW')
        self.add(mb.Particle(name='H', pos=[0.9626, 0.8420, 1.2673]), 'HW1')
        self.add(mb.Particle(name='H', pos=[0.9626, 0.6787, 1.2673]), 'HW2')
        self.add_bond((self.OW, self.HW1))
        self.add_bond((self.OW, self.HW2))

if __name__ == '__main__':
    m = H2O()
    m.visualize()

