from __future__ import division

from numpy import pi

import mbuild as mb

from mbuild.lib.atoms import H
from mbuild.lib.surfaces import Betacristobalite
from mbuild.examples.pmpc.brush import Brush


class PMPCLayer(mb.Monolayer):
    """Create a layer of grafted pMPC brushes on a beta-cristobalite surface."""
    def __init__(self, mask, tile_x=1, tile_y=1, chain_length=4, alpha=pi/4):
        surface = Betacristobalite()
        brush = Brush(chain_length=chain_length, alpha=alpha)
        hydrogen = H()
        super(PMPCLayer, self).__init__(surface, brush, hydrogen,
                                        mask, tile_x, tile_y)


def main():
    mask = mb.grid_mask_2d(2, 1)
    pmpc_system = PMPCLayer(mask=mask, chain_length=5, alpha=pi/4, tile_x=1, tile_y=2)
    return pmpc_system


if __name__ == "__main__":
    pmpc_layer = main()
    #pmpc_layer.save(filename='brush_layer.mol2')
    pmpc_layer.visualize()
    # pmpc_layer.view_hierarchy()
