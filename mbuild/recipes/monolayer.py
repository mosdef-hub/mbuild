from copy import deepcopy
from warnings import warn

import numpy as np

import mbuild as mb

__all__ = ['Monolayer']


class Monolayer(mb.Compound):
    """A general monolayer recipe.

    Parameters
    ----------
    surface : mb.Compound
        Surface on which the monolayer will be built.
    chains : list of mb.Compounds
        The chains to be replicated and attached to the surface.
    fractions : list of floats
        The fractions of the pattern to be allocated to each chain.
    backfill : list of mb.Compound, optional, default=None
        If there are fewer chains than there are ports on the surface,
        copies of `backfill` will be used to fill the remaining ports.
    pattern : mb.Pattern, optional, default=mb.Random2DPattern
        An array of planar binding locations. If not provided, the entire
        surface will be filled with `chain`.
    tile_x : int, optional, default=1
        Number of times to replicate substrate in x-direction.
    tile_y : int, optional, default=1
        Number of times to replicate substrate in y-direction.

    """

    def __init__(self, surface, chains, fractions=None, backfill=None, pattern=None,
                 tile_x=1, tile_y=1, **kwargs):
        super(Monolayer, self).__init__()

        # Replicate the surface.
        tiled_compound = mb.TiledCompound(surface, n_tiles=(tile_x, tile_y, 1))
        self.add(tiled_compound, label='tiled_surface')

        if pattern is None:  # Fill the surface.
            pattern = mb.Random2DPattern(len(tiled_compound.referenced_ports()))

        if isinstance(chains, mb.Compound):
            chains = [chains]
        
        if fractions:
            fractions = list(fractions)
            if len(chains) != len(fractions):
                raise ValueError("Number of fractions does not match the number"
                                 " of chain types provided")

            n_chains = len(pattern.points)

            # Attach chains of each type to binding sites based on
            # respective fractions.
            for chain, fraction in zip(chains[:-1], fractions[:-1]):

                # Create sub-pattern for this chain type
                subpattern = deepcopy(pattern)
                n_points = int(round(fraction * n_chains))
                warn("\n Adding {} of chain {}".format(n_points, chain))
                pick = np.random.choice(subpattern.points.shape[0], n_points,
                                        replace=False)
                points = subpattern.points[pick]
                subpattern.points = points

                # Remove now-occupied points from overall pattern
                pattern.points = np.array([point for point in pattern.points.tolist()
                                           if point not in subpattern.points.tolist()])

                # Attach chains to the surface
                attached_chains, _ = subpattern.apply_to_compound(
                    guest=chain, host=self['tiled_surface'], backfill=None, **kwargs)
                self.add(attached_chains)

        else:
            warn("\n No fractions provided. Assuming a single chain type.")

        # Attach final chain type. Remaining sites get a backfill.
        warn("\n Adding {} of chain {}".format(len(pattern), chains[-1]))
        attached_chains, backfills = pattern.apply_to_compound(guest=chains[-1],
                         host=self['tiled_surface'], backfill=backfill, **kwargs)
        self.add(attached_chains)
        self.add(backfills)
