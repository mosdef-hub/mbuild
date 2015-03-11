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
        self.add(tc, 'tiled_surface')

        brush = Brush(chain_length=chain_length, alpha=alpha)
        hydrogen = H()

        mb.apply_mask(self.tiled_surface, brush, mask, backfill=hydrogen)

def main():
    mask = mb.random_mask_2d(1)
    return PMPCLayer(mask=mask, chain_length=3, alpha=pi/4, tile_x=1, tile_y=1)

if __name__ == "__main__":
    pmpc_layer = main()
    pmpc_layer.visualize(show_ports=True)
    pmpc_layer = pmpc_layer.save(filename='brush_layer.hoomdxml',
                                 chain_types=[Betacristobalite, Brush],
                                 residue_types=[MPC])
