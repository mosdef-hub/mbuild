from __future__ import division

from numpy import pi

import mbuild as mb

from mbuild.components.atoms.H import H
from mbuild.components.surfaces.betacristobalite import Betacristobalite
from mbuild.examples.pmpc.brush import Brush
from mbuild.examples.pmpc.mpc import MPC


class PMPCLayer(mb.Compound):
    """Create a layer of grafted pMPC brushes on a beta-cristobalite surface."""
    def __init__(self, mask, tile_x=1, tile_y=1, chain_length=4, alpha=pi/4):
        super(PMPCLayer, self).__init__()

        surface = Betacristobalite()
        tc = mb.TiledCompound(surface, (tile_x, tile_y, 1), kind="tiled_surface")
        box = tc.boundingbox()
        mb.translate(tc, -box.mins)
        self.add(tc, 'tiled_surface')

        brush = Brush(chain_length=chain_length, alpha=alpha)
        hydrogen = H()

        mb.apply_mask(self.tiled_surface, brush, mask, backfill=hydrogen)

def main():
    mask = mb.grid_mask_2d(5, 5)
    #mask = mb.random_mask_2d(1)
    return PMPCLayer(mask=mask, chain_length=5, alpha=pi/4, tile_x=2, tile_y=1)
    #return PMPCLayer(mask=mask, chain_length=1, alpha=pi/4, tile_x=1, tile_y=1)

if __name__ == "__main__":
    from copy import deepcopy

    import numpy as np

    pmpc_layer = main()
    # box = pmpc_layer.tiled_surface.boundingbox()
    #
    # top_layer = deepcopy(pmpc_layer)
    # mb.rotate_around_x(top_layer, np.pi)
    # mb.translate(top_layer, [0, box.lengths[1], box.lengths[2] * 2.2])
    #
    #
    #
    # out = mb.Compound()
    # out.add(pmpc_layer, 'bot')
    # out.add(top_layer, 'top')
    # out.visualize(show_ports=True)
    pmpc_layer.save(filename='brush_layer.pdb')
    #                             chain_types=[Betacristobalite, Brush],
    #                             residue_types=[MPC])
