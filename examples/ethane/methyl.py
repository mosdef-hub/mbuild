__author__ = 'sallai'
import pdb

from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.mol2file import load_mol2
from mbuild.coordinate_transform import *

class Methyl(Compound):
    """ """
    def __init__(self):
        Compound.__init__(self)
        load_mol2('methyl.mol2', component=self)
        carbon = self.labels['C.3_1']

        transform(self, Translation(-carbon.pos))

        up = Port()
        self.add(up, 'up')
        transform(up, Translation(np.array([0,-0.7,0])))

        down = Port()
        self.add(down, 'down')
        transform(down, RotationAroundZ(np.pi))
        transform(down, Translation(np.array([0,-0.7,0])))

if __name__ == '__main__':
    methyl = Methyl()

    from mbuild.plot import Plot
    Plot(methyl, verbose=True, atoms=True, bonds=True, angles=False, dihedrals=False).show()


