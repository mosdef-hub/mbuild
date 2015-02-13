from __future__ import division

from numpy import pi

from mbuild.compound import Compound
from mbuild.components.surfaces.betacristobalite import Betacristobalite
from mbuild.examples.pmpc_brush_layer.mpc_monomer import MpcMonomer
from mbuild.examples.pmpc_brush_layer.brush import Brush
from mbuild.tools.mask import apply_mask, random_mask_2d
from mbuild.tools.tiled_compound import TiledCompound


class BrushLayer(Compound):
    """Create a layer of grafted pMPC brushes on a beta-cristobalite surface."""
    def __init__(self, tile_x=1, tile_y=1, chain_length=4, alpha=pi/4, mask=None):
        super(BrushLayer, self).__init__()

        surface = Betacristobalite()

        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")
        self.add(tc, 'tiled_surface')
        self.periodicity = tc.periodicity

        brush_proto = Brush(chain_length=chain_length, alpha=alpha)

        if mask is None:
            mask == random_mask_2d(5)
        apply_mask(self.tiled_surface, brush_proto, mask)


def main():
    mask = random_mask_2d(1)
    brush_layer = BrushLayer(chain_length=3, alpha=pi/4, mask=mask, tile_x=1, tile_y=1)

    brush_layer = brush_layer.to_trajectory(chain_types=[Betacristobalite, Brush],
                                            residue_types=[MpcMonomer])
    brush_layer.topology.find_forcefield_terms()
    brush_layer.save(filename='brush_layer.top')

if __name__ == "__main__":
    main()
