from __future__ import division

from mbuild.atom import Atom
from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.tiled_compound import TiledCompound
from mbuild.plugins.mask import apply_mask

from mbuild.examples.alkane_monolayers.alkylsilane import AlkylSilane
from mbuild.examples.alkane_monolayers.surface import Surface


class AlkaneMonolayer(Compound):
    """An akylsilane monolayer on beta-cristobalite. """

    def __init__(self, tile_x=1, tile_y=1, chain_length=4, mask=None):
        """Create an alkylsilane monolayer on beta-cristobalite.

        Args:
            tile_x (int):
            tile_y (int):
            chain_length (int):
            mask (np.ndarray):
        """
        super(AlkaneMonolayer, self).__init__()

        surface = Surface()
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")

        self.add(tc, 'tiled_surface')
        self.periodicity = tc.periodicity

        alkylsilane = AlkylSilane(chain_length)
        # TODO: Library of all elements with ports attached.
        hydrogen = Compound('H')
        hydrogen.add(Atom('H'), 'H')
        hydrogen.add(Port(anchor=hydrogen.H), 'down')
        apply_mask(self.tiled_surface, alkylsilane, mask,
                   guest_port_name='down', backfill=hydrogen)


if __name__ == "__main__":
    from mbuild.plugins.mask import grid_mask_2d
    # mask = random_mask_2d(4)
    mask = grid_mask_2d(3, 3)

    monolayer = AlkaneMonolayer(chain_length=6, mask=mask)
    monolayer = monolayer.to_trajectory()

    import cProfile, pstats, StringIO
    pr = cProfile.Profile()
    pr.enable()

    monolayer.top.find_forcefield_terms()

    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()

    print monolayer.n_atoms
    print monolayer.top.n_ff_bonds
    print monolayer.top.n_ff_angles
    print monolayer.top.n_ff_dihedrals

    #from mbuild.plot import Plot
    #Plot(monolayer, bonds=True, angles=False, dihedrals=False,
    #     periodic_bonds=False).show()

