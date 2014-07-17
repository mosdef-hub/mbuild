from copy import deepcopy

from mbuild.compound import Compound
from mbuild.coordinate_transform import Translation
from mbuild.port import Port
from mbuild.prototype import Prototype
import numpy as np

__author__ = 'sallai'

class TiledCompound(Compound):
    """ """

    def __init__(self, tile, n_x=1, n_y=1, n_z=1, kind=None, ctx={}, label=None):
        assert(isinstance(tile, Compound))
        super(TiledCompound, self).__init__(ctx=ctx)

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

                    #video
                    # new_tile.transform(Translation((i*tile.periodicity[0], j*tile.periodicity[1], 0)))
                    new_tile.transform(Translation((i*(tile.periodicity[0]+8), j*(tile.periodicity[1]+8), 0)))

                    self.add(new_tile,label=label + "_" + str(i)+"_"+str(j), inherit_periodicity=False)
                    for port in new_tile.parts:
                        if isinstance(port, Port):
                            self.add(port, containment=False)

        self.periodicity = np.array([ tile.periodicity[0]*n_x, tile.periodicity[1]*n_y, tile.periodicity[2]*n_z ])


if __name__ == "__main__":

    from mbuild.components.surface import Surface
    surface = Surface()

    Prototype('o-si', color='grey')

    tc = TiledCompound(surface, 3, 4, 1, kind="tiled_surface")

    from mbuild.plot import Plot
    Plot(tc, bonds=True, verbose=True).show()
