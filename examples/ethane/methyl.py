__author__ = 'sallai'
import os, sys
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.file_formats.mol2file import load_mol2
from mbuild.coordinate_transform import *

class Methyl(Compound):
    """ """
    def __init__(self):
        Compound.__init__(self)

        # Look for data file in same directory as this python module.
        # current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        # new_path = os.path.join(current_dir, 'methyl.pdb')
        # load_mol2(new_path, part=self)
        self.append_from_file('methyl.pdb', relative_to_module=__name__)

        # transform(self, Translation(-carbon.pos))
        translate(self, -self.C[0])

        self.add(Port(anchor=self.C[0]), 'up')
        rotate_around_z(self.up, np.pi)
        translate(self.up, np.array([0,-0.7,0]))

        self.add(Port(anchor=self.C[0]), 'down')
        translate(self.down, np.array([0,-0.7,0]))

if __name__ == '__main__':
    methyl = Methyl()

    from mbuild.plot import Plot
    Plot(methyl, verbose=True, atoms=True, bonds=True, angles=False, dihedrals=False).show()


