from __future__ import division

from numpy import pi

from mbuild.compound import Compound
from mbuild.examples.pmpc.surface import Surface
from mbuild.examples.pmpc.brush import Brush
from mbuild.plugins.mask import grid_mask_2d, apply_mask
from mbuild.tiled_compound import TiledCompound


class BrushLayer(Compound):
    """ """
    def __init__(self, tile_x=1, tile_y=1, chain_length=4, alpha=pi/4, mask=None):
        """ """
        super(BrushLayer, self).__init__()

        surface = Surface()
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")
        self.add(tc, 'tiled_surface')
        self.periodicity = tc.periodicity

        brush_proto = Brush(chain_length=chain_length, alpha=alpha)

        apply_mask(self.tiled_surface, brush_proto, mask)

if __name__ == "__main__":
    # mask = random_mask_2d(4)
    mask = grid_mask_2d(3,3)
    # print mask

    m = BrushLayer(chain_length=5, alpha=pi/4, mask=mask, tile_x=2, tile_y=2)

    from mbuild.plot import Plot
    Plot(m, bonds=True, angles=False, dihedrals=False, periodic_bonds=True).show()

    # from mbuild.treeview import TreeView
    # tv = TreeView(m)
    # tv.show()

