from __future__ import division

from numpy import pi

import mbuild as mb

from mbuild.lib.atoms import H
from mbuild.lib.surfaces import Betacristobalite
from mbuild.examples.pmpc.brush import Brush


class PMPCLayer(mb.Monolayer):
    """Create a layer of grafted pMPC brushes on a beta-cristobalite surface."""
    def __init__(self, pattern, tile_x=1, tile_y=1, chain_length=4, alpha=pi / 4):
        surface = Betacristobalite()
        brush = Brush(chain_length=chain_length, alpha=alpha)
        hydrogen = H()
        super(PMPCLayer, self).__init__(surface, brush, backfill=hydrogen,
                                        pattern=pattern, tile_x=tile_x,
                                        tile_y=tile_y)


def main():
    pattern = mb.Random2DPattern(2)
    pmpc_system = PMPCLayer(pattern=pattern, chain_length=3, alpha=pi / 4, tile_x=1, tile_y=1)
    return pmpc_system


if __name__ == "__main__":
    pmpc_layer = main()
    print(pmpc_layer)
