from mbuild.polymer import Polymer
from examples.ethane.methyl import Methyl

__author__ = 'sallai'
from copy import deepcopy
import operator
import pdb

from numpy import pi
from mbuild.bond import Bond

from mbuild.coordinate_transform import *
from mbuild.compound import Compound
from mbuild.port import Port

from ch2 import Ch2

class Alkane(Compound):

    def __init__(self, n=3):
        if n < 3:
            raise Exception('n must be 1 or more')
        Compound.__init__(self)

        self.add(Methyl(), "methyl1")
        self.add( Polymer(Ch2(), n=n-2, port_labels=("up","down")), "chain" )
        self.add(Methyl(), "methyl2")

        equivalence_transform(self.chain, self.chain.up, self.methyl1.down)
        equivalence_transform(self.methyl2, self.methyl2.up, self.chain.down)

if __name__ == "__main__":
    m = Alkane(n=4)
    from mbuild.plot import Plot
    Plot(m).show()
