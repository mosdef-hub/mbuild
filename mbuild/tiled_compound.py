from copy import deepcopy

import numpy as np

from mbuild.bond import Bond
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import translate
from mbuild.periodic_kdtree import PeriodicCKDTree

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

        n_x, n_y, n_z = n_tiles
        if not n_x > 0 and n_y > 0 and n_z > 0:
            raise ValueError("Number of tiles must be positive.")

        # Check that the tile is periodic in the requested dimensions.
        if ((n_x != 1 and tile.periodicity[0] == 0) or
                (n_y != 1 and tile.periodicity[1] == 0) or
                (n_z != 1 and tile.periodicity[2] == 0)):
            raise ValueError("Tile not periodic in at least one of the "
                             "specified dimensions.")

        self.tile = tile
        self.n_x = n_x
        self.n_y = n_y
        self.n_z = n_z
        if kind is None:
            kind = tile.kind
        self.kind = kind

        self._replicate_tiles()
        self._stitch_bonds()
        self.periodicity = np.array([tile.periodicity[0] * n_x,
                                     tile.periodicity[1] * n_y,
                                     tile.periodicity[2] * n_z])

    def _replicate_tiles(self):
        """Replicate and place periodic tiles. """
        for i in range(self.n_x):
            for j in range(self.n_y):
                for k in range(self.n_z):
                    new_tile = deepcopy(self.tile)
                    translate(new_tile,
                              np.array([i * self.tile.periodicity[0],
                                        j * self.tile.periodicity[1],
                                        k * self.tile.periodicity[2]]))
                    tile_label = "{0}_{1}_{2}_{3}".format(self.kind, i, j, k)
                    self.add(new_tile, label=tile_label)

                    # Hoist ports.
                    for port in new_tile.parts:
                        if isinstance(port, Port):
                            self.add(port, containment=False)

    def _stitch_bonds(self):
        """Stitch bonds across periodic boundaries. """
        if not np.all(self.periodicity == 0):
            # For every tile, we assign temporary ID's to atoms.
            for child in self.parts:
                # Not using isinstance because we want to ignore inheritance.
                if type(child) == type(self.tile):
                    for idx, atom in enumerate(child.yield_atoms()):
                        atom.index = idx

            # Build a kdtree of all atoms.
            atom_kdtree = PeriodicCKDTree(self.xyz, bounds=self.periodicity)

            # Cutoff for long bonds is half the shortest periodic distance.
            bond_dist_thres = min(self.tile.periodicity[self.tile.periodicity > 0]) / 2

            # Update connectivity.
            bonds_to_remove = set()
            bonds_to_add = set()
            all_atoms = np.asarray(self.atoms)
            for bond in self.yield_bonds():
                if bond.distance() > bond_dist_thres:
                    # Find new pair for atom1.
                    _, idxs = atom_kdtree.query(bond.atom1.pos, k=10)
                    neighbors = all_atoms[idxs]

                    for atom in neighbors:
                        if atom.index == bond.atom2.index:
                            atom2_image = atom
                            break
                    else:
                        raise RuntimeError('Unable to find matching atom image'
                                           'while stitching bonds.')

                    # Find new pair for atom2.
                    _, idxs = atom_kdtree.query(bond.atom2.pos, k=10)
                    neighbors = all_atoms[idxs]

                    for atom in neighbors:
                        if atom.index == bond.atom1.index:
                            atom1_image = atom
                            break
                    else:
                        raise RuntimeError('Unable to find matching atom image'
                                           'while stitching bonds.')

                    # Mark old bond for removal and add new ones.
                    bonds_to_remove.add(bond)
                    bond1 = Bond(bond.atom1, atom2_image)
                    bonds_to_add.add(bond1)
                    bond2 = Bond(atom1_image, bond.atom2)
                    bonds_to_add.add(bond2)

            # Remove all marked bonds.
            self.remove(bonds_to_remove)
            self.add(bonds_to_add)

            # Remove the temporary index field from all atoms.
            for child in self.parts:
                if child.__class__ != self.tile.__class__:
                    continue
                for atom in child.yield_atoms():
                    del atom.index
