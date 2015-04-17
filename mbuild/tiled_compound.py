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
    def __init__(self, tile, n_tiles=None, kind=None):
        super(TiledCompound, self).__init__()

        if not n_tiles:
            n_tiles = (1, 1, 1)
        n_x, n_y, n_z = n_tiles
        assert n_x > 0 and n_y > 0 and n_z > 0, "Number of tiles must be positive."

        # Check that the tile is periodic in the requested dimensions.
        if n_x != 1:
            assert tile.periodicity[0] != 0, "Tile not periodic in x dimension."
        if n_y != 1:
            assert tile.periodicity[1] != 0, "Tile not periodic in y dimension."
        if n_z != 1:
            assert tile.periodicity[2] != 0, "Tile not periodic in z dimension."

        self.tile = tile
        self.n_x = n_x
        self.n_y = n_y
        self.n_z = n_z
        if kind is None:
            kind = tile.kind
        self.kind = kind

        self.replicate_tiles()
        self.stitch_bonds()
        self.periodicity = np.array([tile.periodicity[0] * n_x,
                                     tile.periodicity[1] * n_y,
                                     tile.periodicity[2] * n_z])
        # TODO: Look into returning a new compound instead of adding the current one.

    def replicate_tiles(self):
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

    def stitch_bonds(self):
        """Stitch bonds across periodic boundaries. """
        if not np.all(self.periodicity == 0):
            # For every tile, we assign temporary ID's to atoms.
            for child in self.parts:
                if child.__class__ != self.tile.__class__:
                    continue
                for idx, atom in enumerate(child.yield_atoms()):
                    atom.uid = idx

            # Build a kdtree of all atoms.
            atom_list = np.array([atom for atom in self.yield_atoms()])
            atom_kdtree = PeriodicCKDTree([atom.pos for atom in atom_list],
                                          bounds=self.periodicity)

            # Cutoff for long bonds is half the shortest periodic distance.
            bond_dist_thres = min(self.tile.periodicity[self.tile.periodicity > 0]) / 2

            # Update connectivity.
            bonds_to_remove = set()
            for bond in self.yield_bonds():
                if bond.distance(self.periodicity) > bond_dist_thres:
                    # Find new pair for atom1.
                    _, idxs = atom_kdtree.query(bond.atom1.pos, k=10)
                    neighbors = atom_list[idxs]

                    atom2_image = None
                    for atom in neighbors:
                        if atom.uid == bond.atom2.uid:
                            atom2_image = atom
                            break
                    assert atom2_image is not None

                    # Find new pair for atom2.
                    _, idxs = atom_kdtree.query(bond.atom2.pos, k=10)
                    neighbors = atom_list[idxs]

                    atom1_image = None
                    for a in neighbors:
                        if a.uid == bond.atom1.uid:
                            atom1_image = a
                            break
                    assert atom2_image is not None

                    # Mark old bond for removal and add new ones.
                    bonds_to_remove.add(bond)
                    self.add(Bond(bond.atom1, atom2_image))
                    self.add(Bond(atom1_image, bond.atom2))

            # Remove all marked bonds.
            self.remove(bonds_to_remove)

            # Remove the temporary uid field from all atoms.
            for child in self.parts:
                if child.__class__ != self.tile.__class__:
                    continue
                for atom in child.yield_atoms():
                    del atom.uid
