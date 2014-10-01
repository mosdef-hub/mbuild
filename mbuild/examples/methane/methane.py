import numpy as np

from mbuild.compound import Compound
from mbuild.atom import Atom
from mbuild.bond import Bond


class Methane(Compound):
    def __init__(self):
        super(Methane, self).__init__()
        self.add(Atom('C', pos=np.array([0, 0, 0])), 'c')
        self.add(Atom('H', pos=np.array([.15, 0, 0])), 'h1')
        self.add(Atom('H', pos=np.array([0, .15, 0])), 'h2')
        self.add(Atom('H', pos=np.array([-.15, 0, 0])), 'h3')
        self.add(Atom('H', pos=np.array([0, -.15, 0])), 'h4')
        self.add(Bond(self.c, self.h1))
        self.add(Bond(self.c, self.h2))
        self.add(Bond(self.c, self.h3))
        self.add(Bond(self.c, self.h4))


def main():
    methane = Methane()
    methane.visualize()

if __name__ == "__main__":
    main()
