from pybel import Molecule

__author__ = 'sallai'
from mbuild.plot import Plot
from mbuild.prototype import Prototype
from mbuild.compound import *
from mbuild.treeview import *

class Methane(Compound):
    def __init__(self):
        super(Methane, self).__init__()
        self.add(Atom('C', pos=np.array([0, 0, 0])), 'c')
        self.add(Atom('H', pos=np.array([.15, 0, 0])), 'h1')
        self.add(Atom('H', pos=np.array([0, .15, 0])), 'h2')
        self.add(Atom('H', pos=np.array([-.15, 0, 0])), 'h3')
        self.add(Atom('H', pos=np.array([0, -.15, 0])), 'h4')
        self.add(Bond(self.c,self.h1))
        self.add(Bond(self.c,self.h2))
        self.add(Bond(self.c,self.h3))
        self.add(Bond(self.c,self.h4))

if __name__ == "__main__":

    m = Methane()

    Prototype("C", radius=.17, color="yellow")
    # Prototype("C", radius=.17, color="black")
    Prototype("H", radius=.12, color="white")
    Prototype("c-h", color="pink")
    Prototype("h-c-h", color="red")
    Prototype("h-h-h-h", color="orange")

    # TreeView(m).show()

    Plot(m, bonds=True, angles=False, dihedrals=False).show()

