import mbuild as mb

from mbuild.components.surfaces.betacristobalite import Betacristobalite
from mbuild.components.atoms.H import H
from mbuild.examples.alkane_monolayer.alkylsilane import AlkylSilane


class AlkaneMonolayer(mb.Compound):
    """An akylsilane monolayer on beta-cristobalite. """

    def __init__(self, mask, tile_x=1, tile_y=1, chain_length=10):
        """Create an alkylsilane monolayer on beta-cristobalite.

        Args:
            tile_x (int): Number of times to replicate substrate in x-direction.
            tile_y (int): Number of times to replicate substrate in y-direction.
            chain_length (int): Number of carbon atoms per chain.
            mask (np.ndarray): A 2D array of binding locations.
        """
        super(AlkaneMonolayer, self).__init__()

        surface = Betacristobalite()
        # Replicate the surface.
        tc = mb.TiledCompound(surface, n_tiles=(tile_x, tile_y, 1), kind="tiled_surface")
        self.add(tc, 'tiled_surface')

        alkylsilane = AlkylSilane(chain_length)
        hydrogen = H()

        # Attach chains to specified binding sites. Other sites get a hydrogen.
        mb.apply_mask(host=self.tiled_surface, guest=alkylsilane, mask=mask,
                   backfill=hydrogen)


def main():
    from mbuild.mask import grid_mask_2d
    mask = grid_mask_2d(8, 8)  # Evenly spaced, 2D grid of points.
    monolayer = AlkaneMonolayer(chain_length=10, mask=mask)

    monolayer.save(filename='c10-n64.pdb', show_ports=True)

if __name__ == "__main__":
    main()

