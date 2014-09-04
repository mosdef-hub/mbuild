import numpy as np
import os, sys

from mbuild.coordinate_transform import *
from mbuild.compound import Compound
from mbuild.file_formats.mol2file import load_mol2
from mbuild.port import Port
from mbuild.tiled_compound import TiledCompound


class Surface(Compound):
    """
    """
    def __init__(self):
        super(Surface, self).__init__()

        # Look for mol2 file in same directory as this file.
        self.append_from_file('beta-cristobalite.mpdb', relative_to_module=__name__)

        self.periodicity = np.array([47.689, 41.3, 0.0])

        cnt = 0
        for atom in self.atoms():
            if atom.kind == 'O' and atom.pos[2] > 10:
                cnt += 1
                port = Port(anchor=atom)
                rotate_around_x(port, pi/2)
                translate(port, atom + np.array([0, 0, 1]))
                self.add(port, 'port_{}'.format(cnt))

if __name__ == "__main__":
    s = Surface()
    m = TiledCompound(s, n_x=2, n_y=1, n_z=1, kind="tiled")

    # print s.periodicity
    # print m.periodicity

    from mbuild.plot import Plot
    Plot(m, bonds=True, verbose=True, periodic_bonds=True).show()
