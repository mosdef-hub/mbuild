from copy import deepcopy

from mbuild.compound import Compound
from mbuild.bond import Bond
from mbuild.coordinate_transform import translate
from mbuild.periodic_kdtree import PeriodicCKDTree

import numpy as np
from mbuild.port import Port

__author__ = 'sallai'

class TiledCompound(Compound):
    """ """

    def __init__(self, tile, n_x=1, n_y=1, n_z=1, kind=None, label=None):

        assert(isinstance(tile, Compound))
        super(TiledCompound, self).__init__()

        assert n_x>0 and n_y>0 and n_z>0, "number or tiles must be positive"

        # check that the tile is periodic in the requested dimensions
        if n_x != 1:
            assert tile.periodicity[0] != 0, "tile is not periodic along the x dimension"
        if n_y != 1:
            assert tile.periodicity[1] != 0, "tile is not periodic along the y dimension"
        if n_z != 1:
            assert tile.periodicity[2] != 0, "tile is not periodic along the z dimension"

        if label is None:
            label = tile.kind

        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    new_tile = deepcopy(tile)

                    translate(new_tile, np.array([i*tile.periodicity[0], j*tile.periodicity[1],  k*tile.periodicity[2]]))

                    self.add(new_tile, label=label + "_" + str(i)+"_"+str(j), inherit_periodicity=False)

                    # hoist ports
                    for port in new_tile.parts:
                        if isinstance(port, Port):
                            self.add(port, containment=False)

        self.periodicity = np.array([ tile.periodicity[0]*n_x, tile.periodicity[1]*n_y, tile.periodicity[2]*n_z ])

        # stitch bonds across periodic boundaries
        if(any(self.periodicity > 0)):

            # for every tile, we assign uids to atoms
            for child in self.parts:
                if child.__class__ != tile.__class__:
                    continue
                for idx, atom in enumerate(child.atoms()):
                    atom.uid = idx

            # build a kdtree of all atoms
            atomsList = np.array([atom for atom in self.atoms()])
            atomKdtree = PeriodicCKDTree([atom.pos for atom in atomsList], bounds=self.periodicity)

            # for every long bond
            bond_dist_thres = min(tile.periodicity[tile.periodicity > 0]) / 2

            print "Stitching bonds..."
            bonds_to_remove = set()
            for bond in self.bonds():
                if bond.distance(self.periodicity) > bond_dist_thres:

                    # find new pair for atom1
                    dists, idxs = atomKdtree.query(bond.atom1.pos, k=10)
                    neighbors = atomsList[idxs]

                    atom2_image = None
                    for a in neighbors:
                        if a.uid == bond.atom2.uid:
                            atom2_image = a
                            break

                    assert(atom2_image is not None)

                    # find new pair for atom2
                    dists, idxs = atomKdtree.query(bond.atom2.pos, k=10,)
                    neighbors = atomsList[idxs]

                    atom1_image = None
                    for a in neighbors:
                        if a.uid == bond.atom1.uid:
                            atom1_image = a
                            break

                    assert(atom2_image is not None)

                    # remove existing bond
                    bonds_to_remove.add(bond)
                    # add bond between atom1 and atom2_image
                    self.add(Bond(bond.atom1, atom2_image))
                    # add bond between atom1_image and atom2
                    self.add(Bond(atom1_image, bond.atom2))

            # batch remove bonds
            self.remove(bonds_to_remove, containment_only=True)

            # for every tile, we clean up uids
            for child in self.parts:
                if child.__class__ != tile.__class__:
                    continue
                for atom in child.atoms():
                    del atom.uid
            print "Stitching bonds done..."

if __name__ == "__main__":

    from examples.pmpc.surface import Surface

    surface = Surface()

    tc = TiledCompound(surface, 2, 3, 1, kind="tiled_surface")

    from mbuild.plot import Plot
    Plot(tc, bonds=True, verbose=True, periodic_bonds=True).show()
