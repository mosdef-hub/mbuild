from __future__ import division

from mbuild.compound import Compound
from mbuild.coordinate_transform import equivalence_transform
from mbuild.examples.alkane.ch2 import Ch2
from mbuild.examples.alkane_monolayers.alkylsilane import AlkylSilane
from mbuild.examples.ethane.methyl import Methyl
from mbuild.tiled_compound import TiledCompound
from mbuild.plugins.mask import apply_mask

from mbuild.examples.alkane_monolayers.surface import Surface


class AlkaneMonolayer(Compound):
    """ """
    def __init__(self, tile_x=1, tile_y=1, chain_length=4, mask=None):
        """ """
        super(AlkaneMonolayer, self).__init__()

        surface = Surface()
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")

        self.add(tc, 'tiled_surface')
        self.periodicity = tc.periodicity

        #alkylsilane = AlkylSilane(chain_length)
        ch3 = Methyl()
        ch2 = Ch2()
        alkylsilane = Compound()
        alkylsilane.add(ch3, 'ch3')
        alkylsilane.add(ch2, 'ch2')
        equivalence_transform(alkylsilane.ch3, alkylsilane.ch3.down, alkylsilane.ch2.up)
        alkylsilane.add(ch2.down, 'down', containment=False)

        import pdb
        #pdb.set_trace()

        apply_mask(self.tiled_surface, alkylsilane, mask, guest_port_name='down')

if __name__ == "__main__":
    from mbuild.plugins.mask import grid_mask_2d
    # mask = random_mask_2d(4)
    mask = grid_mask_2d(3, 3)

    monolayer = AlkaneMonolayer(chain_length=6, mask=mask)

    from mbuild.plot import Plot
    Plot(monolayer, bonds=True, angles=False, dihedrals=False, periodic_bonds=True).show()

    # from mbuild.treeview import TreeView
    # tv = TreeView(m)
    # tv.show()

