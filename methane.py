from mbuild.plot import Plot
from mbuild.prototype import Prototype

__author__ = 'sallai'
from mbuild.compound import *
from mbuild.treeview import *

class Methane(Compound):
    def __init__(self, ctx={}):
        super(Methane, self).__init__(ctx=ctx)
        self.add(Atom('C', pos=(0, 0, 0)), 'c')
        self.add(Atom('H', pos=(1.5, 0, 0)), 'h1')
        self.add(Atom('H', pos=(0, 1.5, 0)), 'h2')
        self.add(Atom('H', pos=(-1.5, 0, 0)), 'h3')
        self.add(Atom('H', pos=(0, -1.5, 0)), 'h4')

        self.add(Bond(self.c,self.h1,"c-h"))
        self.add(Bond(self.c,self.h2,"c-h"))
        self.add(Bond(self.c,self.h3,"c-h"))
        self.add(Bond(self.c,self.h4,"c-h"))

if __name__ == "__main__":


    m = Methane()

    Prototype("C", radius=1.7, color="yellow")
    # Prototype("C", radius=1.7, color="black")
    # Prototype("H", radius=1.2, color="white")
    Prototype("c-h", color="pink")

    # TreeView(m).show()

    Plot(m, bonds=True).show()

