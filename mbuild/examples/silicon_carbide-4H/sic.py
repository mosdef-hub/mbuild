from __future__ import division

from mbuild.coordinate_transform import *
from mbuild.compound import Compound
from mbuild.tiled_compound import TiledCompound
from mbuild.testing.tools import get_fn


class Surface(Compound):
    """ """
    def __init__(self):
        super(Surface, self).__init__()

        self.append_from_file(get_fn('sic.pdb'))
        atoms_to_remove = list()
        for atom in self.atoms():
            if atom.pos[2] < 0.05:
                atoms_to_remove.append(atom)
        self.remove(atoms_to_remove)
        self.periodicity = np.array([2*6.16 + 1.5, 2*3.56+1.5, 10.08+1.5]) / 10

if __name__ == "__main__":
    single_surface = Surface()
    multi_surface = TiledCompound(single_surface, n_x=8, n_y=8, n_z=2, kind="tiled")

    multi_surface.save('sic-tiled.pdb')
    from mbuild.plot import Plot
    Plot(multi_surface, bonds=True, verbose=True, periodic_bonds=True).show()
