import mbuild as mb


class Monolayer(mb.Compound):
    """A general monolayer recipe.

    Parameters
    ----------
    surface : mb.Compound
        Surface on which the monolayer will be built.
    chain : list of mb.Compound
        The chain to be replicated and attached to the surface.
    backfill : list of mb.Compound, optional, default=None
        If there are fewer chains than there are ports on the surface,
        copies of `backfill` will be used to fill the remaining ports.
    mask : np.ndarray, shape=(n, 3), optional, default=None
        An array of planar binding locations. If not provided, the entire
        surface will be filled with `chain`.
    tile_x : int, optional, default=1
        Number of times to replicate substrate in x-direction.
    tile_y : int, optional, default=1
        Number of times to replicate substrate in y-direction.

    TODO
    ----
    * Support for arbitrary number of chain types and relative proportions

    """

    def __init__(self, surface, chain, backfill=None, mask=None, tile_x=1,
                 tile_y=1, **kwargs):
        super(Monolayer, self).__init__()

        # Replicate the surface.
        tc = mb.TiledCompound(surface, n_tiles=(tile_x, tile_y, 1))
        self.add(tc, 'tiled_surface')

        if mask is None:  # Fill the surface.
            mask = mb.random_mask_2d(len(tc.referenced_ports()))

        # Attach chains to specified binding sites. Remaining sites get a backfill.
        chains, backfills = mb.apply_mask(host=self.tiled_surface, guest=chain,
                                          mask=mask, backfill=backfill, **kwargs)
        self.add(chains)
        self.add(backfills)


if __name__ == "__main__":
    from mbuild.lib.surfaces import Betacristobalite
    from mbuild.examples import AlkylSilane
    from mbuild.lib.atoms import H

    mask = mb.grid_mask_2d(8, 8)  # Evenly spaced, 2D grid of points.
    monolayer = Monolayer(surface=Betacristobalite(), chain=AlkylSilane(10),
                          backfill=H(), mask=mask)
    monolayer.visualize()
