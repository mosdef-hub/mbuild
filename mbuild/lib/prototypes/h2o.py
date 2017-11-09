import mbuild as mb


class H2O(mb.Compound):
    """A water molecule. """
    def __init__(self):
        super(H2O, self).__init__()

        self.add(mb.Particle(name='OW', pos=[1.0203, 0.7604, 1.2673]))
        self.add(mb.Particle(name='HW1', pos=[0.9626, 0.8420, 1.2673]))
        self.add(mb.Particle(name='HW2', pos=[0.9626, 0.6787, 1.2673]))
        self.add_bond((self[0], self[1]))
        self.add_bond((self[0], self[2]))
        self.name = 'HOH'

if __name__ == '__main__':
    m = H2O()
    m.save('h2o.mol2')
