__all__ = ['TiledCompound']

import itertools as it

import numpy as np

from mbuild.compound import Compound
from mbuild.exceptions import MBuildError
from mbuild.port import Port
from mbuild.periodic_kdtree import PeriodicCKDTree
from mbuild import clone


class TiledCompound(Compound):
    """Replicates a Compound in any cartesian direction(s).

    Correctly updates connectivity while respecting periodic boundary
    conditions.

    Parameters
    -----------
    tile : mb.Compound
        The Compound to be replicated.
    n_tiles : array-like, shape=(3,), dtype=int, optional, default=(1, 1, 1)
        Number of times to replicate tile in the x, y and z-directions.
    name : str, optional, default=tile.name
        Descriptive string for the compound.

    """
    def __init__(self, tile, n_tiles, name=None):
        super(TiledCompound, self).__init__()

        n_tiles = np.asarray(n_tiles)
        if not np.all(n_tiles > 0):
            raise ValueError('Number of tiles must be positive.')

        # Check that the tile is periodic in the requested dimensions.
        if np.any(np.logical_and(n_tiles != 1, tile.periodicity == 0)):
            raise ValueError('Tile not periodic in at least one of the '
                             'specified dimensions.')

        if name is None:
            name = tile.name + '-'.join(str(d) for d in n_tiles)
        self.name = name
        self.periodicity = np.array(tile.periodicity * n_tiles)

        if all(n_tiles == 1):
            self._add_tile(tile, [(0, 0, 0)])
            self._hoist_ports(tile)
            return  # Don't waste time copying and checking bonds.

        # For every tile, assign temporary ID's to particles which are internal
        # to that tile. E.g., when replicating a tile with 1800 particles, every
        # tile will contain particles with ID's from 0-1799. These ID's are used
        # below to fix bonds crossing periodic boundary conditions where a new
        # tile has been placed.
        for idx, particle in enumerate(tile.particles(include_ports=True)):
            particle.index = idx

        # Replicate and place periodic tiles.
        # -----------------------------------
        for ijk in it.product(range(n_tiles[0]),
                              range(n_tiles[1]),
                              range(n_tiles[2])):
            new_tile = clone(tile)
            new_tile.translate(np.array(ijk * tile.periodicity))
            self._add_tile(new_tile, ijk)
            self._hoist_ports(new_tile)

        # Fix bonds across periodic boundaries.
        # -------------------------------------
        # Cutoff for long bonds is half the shortest periodic distance.
        bond_dist_thres = min(tile.periodicity[tile.periodicity > 0]) / 2

        # Bonds that were periodic in the original tile.
        indices_of_periodic_bonds = set()
        for particle1, particle2 in tile.bonds():
            if np.linalg.norm(particle1.pos - particle2.pos) > bond_dist_thres:
                indices_of_periodic_bonds.add((particle1.index,
                                                        particle2.index))

        # Build a periodic kdtree of all particle positions.
        self.particle_kdtree = PeriodicCKDTree(data=self.xyz, bounds=self.periodicity)
        all_particles = np.asarray(list(self.particles(include_ports=False)))

        # Store bonds to remove/add since we'll be iterating over all bonds.
        bonds_to_remove = set()
        bonds_to_add = set()
        for particle1, particle2 in self.bonds():
            if (particle1.index, particle2.index) not in indices_of_periodic_bonds \
                    and (particle2.index, particle1.index) not in indices_of_periodic_bonds:
                continue

            if self.min_periodic_distance(particle1.pos, particle2.pos) > bond_dist_thres:
                bonds_to_remove.add((particle1, particle2))
                image2 = self._find_particle_image(particle1, particle2,
                                                   all_particles)
                image1 = self._find_particle_image(particle2, particle1,
                                                   all_particles)

                if (image2, particle1) not in bonds_to_add:
                    bonds_to_add.add((particle1, image2))
                if (image1, particle2) not in bonds_to_add:
                    bonds_to_add.add((particle2, image1))

        for bond in bonds_to_remove:
            self.remove_bond(bond)

        for bond in bonds_to_add:
            self.add_bond(bond)

        # Clean up temporary data.
        for particle in self._particles(include_ports=True):
            particle.index = None
        del self.particle_kdtree

    def _add_tile(self, new_tile, ijk):
        """Add a tile with a label indicating its tiling position. """
        tile_label = "{0}_{1}".format(self.name, '-'.join(str(d) for d in ijk))
        self.add(new_tile, label=tile_label, inherit_periodicity=False)

    def _hoist_ports(self, new_tile):
        """Add labels for all the ports to the parent (TiledCompound). """
        for port in new_tile.children:
            if isinstance(port, Port):
                self.add(port, containment=False)

    def _find_particle_image(self, query, match, all_particles):
        """Find particle with the same index as match in a neighboring tile. """
        _, idxs = self.particle_kdtree.query(query.pos, k=10)

        neighbors = all_particles[idxs]

        for particle in neighbors:
            if particle.index == match.index:
                return particle
        raise MBuildError('Unable to find matching particle image while'
                          ' stitching bonds.')
