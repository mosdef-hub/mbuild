from __future__ import division

from numpy import pi

import mbuild as mb

from mbuild.components.atoms.H import H
from mbuild.components.surfaces.betacristobalite import Betacristobalite
from mbuild.examples.pmpc.brush import Brush


class PMPCLayer(mb.Compound):
    """Create a layer of grafted pMPC brushes on a beta-cristobalite surface."""
    def __init__(self, mask, tile_x=1, tile_y=1, chain_length=4, alpha=pi/4):
        super(PMPCLayer, self).__init__()

        surface = Betacristobalite()
        tc = mb.TiledCompound(surface, (tile_x, tile_y, 1))
        box = tc.boundingbox
        mb.translate(tc, -box.mins)
        self.add(tc, 'tiled_surface')

        brush = Brush(chain_length=chain_length, alpha=alpha)
        hydrogen = H()

        mb.apply_mask(self.tiled_surface, brush, mask, backfill=hydrogen)


def main():
    mask = mb.grid_mask_2d(5, 5)
    pmpc_layer = PMPCLayer(mask=mask, chain_length=5, alpha=pi/4, tile_x=2, tile_y=2)
    return pmpc_layer


if __name__ == "__main__":

    pmpc_layer = main()
    pmpc_layer.save(filename='brush_layer.mol2')
    #                             chain_types=[Betacristobalite, Brush],
    #                             residue_types=[MPC])
