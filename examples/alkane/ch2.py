__author__ = 'sallai'

from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.file_formats.mol2file import load_mol2
from mbuild.coordinate_transform import *

class Ch2(Compound):
    """ """
    def __init__(self):
        Compound.__init__(self)

        # Look for data file in same directory as this python module.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        new_path = os.path.join(current_dir, 'ch2.mol2')
        load_mol2(new_path, component=self)

        self.add(Port(anchor=self.C_1), 'up')
        translate(self.up, np.array([0,0.7,0]))

        self.add(Port(anchor=self.C_1), 'down')
        translate(self.down, np.array([0,-0.7,0]))

if __name__ == '__main__':
    m = Ch2()

    from mbuild.plot import Plot
    Plot(m, verbose=True, atoms=True, bonds=True, angles=False, dihedrals=False).show()


