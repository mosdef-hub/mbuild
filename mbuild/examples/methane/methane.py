from mbuild.compound import Compound
from mbuild.atom import Atom
from mbuild.bond import Bond


class Methane(Compound):
    def __init__(self):
        super(Methane, self).__init__()
        carbon = Atom(kind='C')
        self.add(carbon)

        hydrogen = Atom(kind='H', pos=[0.15, 0, 0])
        self.add(hydrogen, label='hc[$]')

        ch_bond = Bond(self.atoms[0], self.hc[0])
        self.add(ch_bond)

        self.add(Atom(kind='H', pos=[0, 0.15, 0]), 'hc[$]')
        self.add(Bond(self.atoms[0], self.hc[1]))
        self.add(Atom(kind='H', pos=[-0.15, 0, 0]), 'hc[$]')
        self.add(Bond(self.atoms[0], self.hc[2]))
        self.add(Atom(kind='H', pos=[0, -0.15, 0]), 'hc[$]')
        self.add(Bond(self.atoms[0], self.hc[3]))

def main():
    methane = Methane()
    methane.visualize()

if __name__ == "__main__":
    main()
