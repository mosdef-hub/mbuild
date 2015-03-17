import mbuild as mb


class H2O(mb.Compound):
    """A water molecule. """
    def __init__(self):
        super(H2O, self).__init__(self)

        self.add(mb.Atom(name='O', pos=[0, 0, 0]), 'OW')
        self.add(mb.Atom(name='H', pos=[0, 0.02, 0]), 'HW1')
        self.add(mb.Atom(name='H', pos=[-0.02, -0.01, 0]), 'HW2')
        self.add(mb.Bond(self.OW, self.HW1))
        self.add(mb.Bond(self.OW, self.HW2))

if __name__ == '__main__':
    m = H2O()
    m.visualize()

