import mbuild as mb


class Methane(mb.Compound):
    def __init__(self):
        super(Methane, self).__init__()
        carbon = mb.Atom(name='C')
        self.add(carbon, label='C')

        hydrogen = mb.Atom(name='H', pos=[0.1, 0, -0.07])
        self.add(hydrogen, label='HC[$]')

        ch_bond = mb.Bond(self.atoms[0], self.HC[0])
        self.add(ch_bond)

        self.add(mb.Atom(name='H', pos=[-0.1, 0, -0.07]), label='HC[$]')
        self.add(mb.Atom(name='H', pos=[0, 0.1, 0.07]), label='HC[$]')
        self.add(mb.Atom(name='H', pos=[0, -0.1, 0.07]), label='HC[$]')

        self.add(mb.Bond(self.atoms[0], self.HC[1]))
        self.add(mb.Bond(self.atoms[0], self.HC[2]))
        self.add(mb.Bond(self.atoms[0], self.HC[3]))


def main():
    methane = Methane()
    return methane

if __name__ == "__main__":
    methane = main()
    methane.visualize()
