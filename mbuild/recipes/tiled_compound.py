import itertools as it

import numpy as np

from mbuild.bond import Bond
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import translate
from mbuild.periodic_kdtree import PeriodicCKDTree
from mbuild import clone
__all__ = ['TiledCompound']


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
    kind : str, optional, default=tile.kind
        Descriptive string for the compound.

    """
    def __init__(self, tile, n_tiles, kind=None):
        super(TiledCompound, self).__init__()

        n_tiles = np.asarray(n_tiles)
        if not np.all(n_tiles > 0):
            raise ValueError('Number of tiles must be positive.')

        # Check that the tile is periodic in the requested dimensions.
        if np.any(np.logical_and(n_tiles != 1, tile.periodicity == 0)):
            raise ValueError('Tile not periodic in at least one of the '
                             'specified dimensions.')

        if kind is None:
            kind = tile.kind + '-'.join(str(d) for d in n_tiles)
        self.kind = kind
        self.periodicity = np.array(tile.periodicity * n_tiles)

        if all(n_tiles == 1):
            self._add_tile_and_hoist_ports(tile, [(0, 0, 0)])
            return  # Don't waste time copying and checking bonds.

        # For every tile, assign temporary ID's to atoms which are internal to
        # that tile. E.g., when replicating a tile with 1800 atoms, every tile
        # will contain atoms with ID's from 0-1799. These ID's are used below
        # to fix bonds crossing periodic boundary conditions where a new tile
        # has been placed.
        for idx, atom in enumerate(tile.yield_atoms()):
            atom.index = idx

        # Replicate and place periodic tiles.
        # -----------------------------------
        for ijk in it.product(range(n_tiles[0]),
                              range(n_tiles[1]),
                              range(n_tiles[2])):
            new_tile = clone(tile)
            translate(new_tile, np.array(ijk * tile.periodicity))
            self._add_tile_and_hoist_ports(new_tile, ijk)



        # Fix bonds across periodic boundaries.
        # -------------------------------------
        # Cutoff for long bonds is half the shortest periodic distance.
        bond_dist_thres = min(tile.periodicity[tile.periodicity > 0]) / 2

        # Bonds that were periodic in the original tile.
        atom_indices_of_periodic_bonds = set()
        for bond in tile.yield_bonds():
            if bond.length() > bond_dist_thres:
                atom_indices_of_periodic_bonds.add((bond.atom1.index, bond.atom2.index))

        # Build a periodic kdtree of all atom positions.
        self.atom_kdtree = PeriodicCKDTree(data=self.xyz, bounds=self.periodicity)
        all_atoms = np.asarray(self.atoms)

        # Store bonds to remove/add since we'll be iterating over all bonds.
        bonds_to_remove = set()
        bonds_to_add = set()
        for bond in self.yield_bonds():
            atom_indices = (bond.atom1.index, bond.atom2.index)
            if atom_indices in atom_indices_of_periodic_bonds:
                if bond.length(self.periodicity) > bond_dist_thres:
                    bonds_to_remove.add(bond)

                    atom2_image = self._find_atom_image(bond.atom1, bond.atom2, all_atoms)
                    new_bond = Bond(bond.atom1, atom2_image)
                    bonds_to_add.add(new_bond)

        self.remove(bonds_to_remove)
        self.add(bonds_to_add)

        # Clean up temporary data.
        for atom in self.yield_atoms():
            atom.index = None
        del self.atom_kdtree

    def _add_tile_and_hoist_ports(self, new_tile, ijk):
        tile_label = "{0}_{1}".format(self.kind, '-'.join(str(d) for d in ijk))
        self.add(new_tile, label=tile_label, inherit_periodicity=False)

        # Hoist ports.
        for port in new_tile.parts:
            if isinstance(port, Port):
                self.add(port, containment=False)

    def _find_atom_image(self, query, match, all_atoms):
        """Find atom with the same index as match in a neighboring tile. """
        _, idxs = self.atom_kdtree.query(query.pos, k=10)
        neighbors = all_atoms[idxs]

        for atom in neighbors:
            if atom.index == match.index:
                return atom
        else:
            raise RuntimeError('Unable to find matching atom image while '
                               'stitching bonds.')