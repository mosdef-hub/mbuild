__author__ = 'sallai'

from methyl import Methyl
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.mol2file import load_mol2
from mbuild.coordinate_transform import *


class Ethane(Compound):
    """ """

    def __init__(self):
        super(Ethane, self).__init__(kind='Ethane')
        methyl1 = Methyl()
        port1 = methyl1.references['up']

        methyl2 = Methyl()
        port2 = methyl2.references['down']

        transform(methyl2, [(port1, port2)])

if __name__ == "__main__":
    ethane = Ethane()

    from mbuild.plot import Plot
    Plot(methyl, verbose=True, atoms=True, bonds=True, angles=False, dihedrals=False).show()
