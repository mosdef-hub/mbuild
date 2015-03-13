import mbuild as mb


class Methane(mb.Compound):
    def __init__(self):
        super(Methane, self).__init__()
        carbon = mb.Atom(name='C')
        self.add(carbon)

        hydrogen = mb.Atom(name='H', pos=[0.15, 0, 0])
        self.add(hydrogen, label='hc[$]')

        ch_bond = mb.Bond(self.atoms[0], self.hc[0])
        self.add(ch_bond)

        self.add(mb.Atom(name='H', pos=[0, 0.15, 0]), 'hc[$]')
        self.add(mb.Bond(self.atoms[0], self.hc[1]))
        self.add(mb.Atom(name='H', pos=[-0.15, 0, 0]), 'hc[$]')
        self.add(mb.Bond(self.atoms[0], self.hc[2]))
        self.add(mb.Atom(name='H', pos=[0, -0.15, 0]), 'hc[$]')
        self.add(mb.Bond(self.atoms[0], self.hc[3]))

def main():
    methane = Methane()
    return methane

if __name__ == "__main__":
    main()
