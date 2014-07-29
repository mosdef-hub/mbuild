from copy import deepcopy

from compound import Compound
from coordinate_transform import Translation
from mbuild import periodic_kdtree
from mbuild.coordinate_transform import translate
from mbuild.periodic_kdtree import PeriodicCKDTree
from port import Port
from prototype import Prototype
import numpy as np

__author__ = 'sallai'

class TiledCompound(Compound):
    """ """

    def __init__(self, tile, n_x=1, n_y=1, n_z=1, kind=None, label=None):
        assert(isinstance(tile, Compound))
        super(TiledCompound, self).__init__()

        assert n_x>0 and n_y>0 and n_z>0, "number or tiles must be positive"

        # check that the tile is periodic in the requested dimensions
        if n_x != 1:
            assert tile.periodicity[0] != 0, "tile is not periodic along the x dimension"
        if n_y != 1:
            assert tile.periodicity[1] != 0, "tile is not periodic along the y dimension"
        if n_z != 1:
            assert tile.periodicity[2] != 0, "tile is not periodic along the z dimension"

        if label is None:
            label = tile.kind

        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    new_tile = deepcopy(tile)

                    translate(new_tile, np.array([i*tile.periodicity[0], j*tile.periodicity[1],  k*tile.periodicity[2]]))

                    self.add(new_tile, label=label + "_" + str(i)+"_"+str(j), inherit_periodicity=False)
                    # for port in new_tile.parts:
                    #     if isinstance(port, Port):
                    #         self.add(port, containment=False)

        self.periodicity = np.array([ tile.periodicity[0]*n_x, tile.periodicity[1]*n_y, tile.periodicity[2]*n_z ])

        #
        #
        # atomKdtree = PeriodicCKDTree([atom.pos for atom in self.atomsList], bounds=self.periodicity)
        #
        # labels = [label for ]
        #
        #
        # for bond in tile.bonds():
        #     if ((tile.periodicity[0] > 0 and np.abs(bond.atom1.pos[0] - bond.atom2.pos[0]) > tile.periodicity[0] / 2.0)
        #         or (tile.periodicity[1] > 0 and np.abs(bond.atom1.pos[1] - bond.atom2.pos[1]) > tile.periodicity[1] / 2.0 )
        #         or (tile.periodicity[2] > 0 and np.abs(bond.atom1.pos[2] - bond.atom2.pos[2]) > tile.periodicity[2] / 2.0)):
        #         # let's find an atom2 that is closer than the current one
        #         atom2_label = atom2



if __name__ == "__main__":

    from mbuild.components.surface import Surface
    surface = Surface()

    Prototype('o-si', color='grey')

    tc = TiledCompound(surface, 3, 4, 1, kind="tiled_surface")

    from mbuild.plot import Plot
    Plot(tc, bonds=True, verbose=True).show()
