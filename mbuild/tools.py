# __author__ = 'sallai'
from __future__ import division
from copy import deepcopy
from atom import Atom

import numpy as np

from compound import Compound
from bond import Bond
from box import Box
from coordinate_transform import equivalence_transform, translate
from periodic_kdtree import PeriodicCKDTree


def vdw_radius(atomic_number):
    return .2


def solvate(host_compound, guest_compound, host_bounds, guest_bounds):
    """Solvate a Compound in a Box.
   
    Typical usage of this function would be to use a pre-equilibrated solvent
    system as the guest which will then be replicated to fill the host's Box.
    All Compounds that have overlapping Atoms with the host are removed.

    Args:
        host_compound (Compound):
        guest_compound (Compound):
        host_bounds (Box):
        guest_bounds (Box):

    """
    assert isinstance(host_bounds, Box)
    assert isinstance(guest_bounds, Box)
    assert np.all(guest_compound.periodicity == 0)

    # TODO: we may want to make sure that the axes of the two boxes line up

    # Replicate the guest so that it fills or overfills the host box.
    host_atom_list = [atom for atom in host_compound.atoms() if atom.kind != 'G']
    host_atom_pos_list = [atom.pos for atom in host_atom_list]
    kdtree = PeriodicCKDTree(host_atom_pos_list)

    num_replicas = np.ceil(host_bounds.lengths / guest_bounds.lengths)
    num_replicas = num_replicas.astype('int')
    for xi in range(0,num_replicas[0]):
        for yi in range(0,num_replicas[0]):
            for zi in range(0,num_replicas[0]):
                guest = deepcopy(guest_compound)
                translate(guest, -guest_bounds.mins + host_bounds.mins + np.array([xi, yi, zi])*guest_bounds.lengths)

                # Remove atoms outside the host's box and anything bonded to them.
                guest_atoms = guest.atom_list_by_kind('*')
                atoms_to_remove = set()
                guest_atom_pos_list = [atom.pos for atom in guest_compound.atoms()]
                atom_indicies = np.where(np.logical_or(np.any(guest_atom_pos_list < host_bounds.mins, axis=1),
                        np.any(guest_atom_pos_list > host_bounds.maxes, axis=1)))[0]
                for ai in atom_indicies:
                    atoms_to_remove.add(guest_atoms[ai])
                    atoms_to_remove.update(guest_atoms[ai].bonded_atoms())
                guest.remove(atoms_to_remove)

                # Remove overlapping atoms and anything bonded to them.
                atoms_to_remove = set()
                for guest_atom in guest.atoms():
                    _, neighbors = kdtree.query(guest_atom.pos, k=10)
                    for host_atom_idx in neighbors:
                        if host_atom_idx < len(host_atom_list):
                            host_atom = host_atom_list[host_atom_idx]
                            if host_compound.min_periodic_distance(host_atom, guest_atom) < (vdw_radius(host_atom) + vdw_radius(guest_atom)):
                                atoms_to_remove.add(guest_atom)
                                atoms_to_remove.update(guest_atom.bonded_atoms())

                guest.remove(atoms_to_remove)
                host_compound.add(guest, "guest_{0}_{1}_{2}".format(xi,yi,zi))


def add_bond(compound, type_A, type_B, dmin, dmax, kind=None):
    """Ai-Bj distance is in [dmin, dmax] => add bond A1xB(Ai,Bj) (symmetric)."""

    for a1 in compound.atom_list_by_kind(type_A):
        nearest = compound.atoms_in_range(a1.pos, dmax)
        for a2 in nearest:
            if (a2.kind==type_B) and (dmin <= compound.min_periodic_distance(a2.pos, a1.pos) <= dmax):
                compound.add(Bond(a1, a2, kind=kind))

def find_all_angles(traj):
    """ """


def add_angle(traj, type_A, type_B, type_C):
    """
    """
    for ab in traj.bonds_by_atom_type(type_A, type_B):
        nearest = traj.neighbor_bonds(ab)
        for bc in nearest:
            bc = traj.order_bond(bc, type_B, type_C)
            if bc is None:
                continue

            temp_ang = None
            if ab[1] == bc[0]:
                temp_ang = (ab[0], ab[1], bc[1])

            if ab[0] == bc[1]:
                temp_ang = (ab[0], ab[1], bc[1])

            if temp_ang is None:
                continue

            # to be improved...
            # handle an edge case where we're adding the same angle multiple times
            if type_A == type_B == type_C:
                if ab[0].index < bc[1].index:
                    continue

            traj.top.add_angle(*temp_ang)

    #
    # def add_dihedral(self, type_A, type_B, type_C, type_D, dihedralKind):
    #     """
    #     """
    #
    #     for abc1 in self.compound.getAnglesByAtomKind(type_A, type_B, type_C):
    #         abc = abc1.cloneWithOrder(type_A, type_B, type_C)
    #
    #         nearest = self.compound.getAnglesInRange(abc.atom2.pos, 5)
    #         for bcd1 in nearest:
    #             if abc1 == bcd1:
    #                 continue
    #
    #             if not bcd1.hasAtomKinds(type_B, type_C, type_D):
    #                 continue
    #
    #             bcd = bcd1.cloneWithOrder(type_B, type_C, type_D)
    #
    #             temp_dhdr = Dihedral.createFromAngles(abc, bcd, kind=dihedralKind)
    #             if temp_dhdr:
    #                 self.compound.add(temp_dhdr)



if __name__ == "__main__":
    print "hello"



