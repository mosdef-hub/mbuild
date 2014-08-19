from __future__ import division
from copy import deepcopy
from atom import Atom
from compound import Compound
from mbuild.box import Box
from mbuild.coordinate_transform import equivalence_transform, translate
from mbuild.periodic_kdtree import PeriodicCKDTree

__author__ = 'sallai'

import numpy as np


def apply_mask(host, guest, mask):
    box = host.boundingbox(excludeG=False)

    mask = mask * box.lengths + box.mins

    print mask

    n_ports = len(host.referenced_ports())
    assert(n_ports >= mask.shape[0])

    port_pos = np.empty((n_ports,3))
    port_list = []
    for pidx, port in enumerate(host.referenced_ports()):
        port_pos[pidx, :] = port.middle.pos
        port_list.append(port)

    for mp in mask:
        closest_point_idx = np.argmin(host.min_periodic_distance(mp, port_pos))
        closest_port = port_list[closest_point_idx]
        brush = deepcopy(guest)
        equivalence_transform(brush, brush.port, closest_port)
        host.add(brush)
        port_pos[closest_point_idx,:] = np.array([np.inf, np.inf, np.inf])


def random_mask_3d(num_sites):
    mask = np.random.random((num_sites, 3))
    return mask


def random_mask_2d(num_sites):
    mask = random_mask_3d(num_sites)
    mask[:, 2] = 0
    return mask


def grid_mask_2d(n, m):
    mask = np.zeros(shape=(n*m, 3), dtype=float)
    for i in range(n):
        for j in range(m):
            mask[i*m + j, 0] = i / n
            mask[i*m + j, 1] = j / m
    return mask


def vdw_radius(atomic_number):
    return 1.5

def solvate(host_compound, guest_compound, host_bounds, guest_bounds):
    assert(isinstance(host_bounds, Box))
    assert(isinstance(guest_bounds, Box))
    assert(np.all(guest_compound.periodicity == 0))

    # we may want to make sure that the axes of the two boxes line up

    # replicate the quest so that it's bigger than the host
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

                # remove atoms outside the host's box, and anything bonded to them (recursively)
                guest_atoms = guest.getAtomListByKind('*')
                atoms_to_remove = set()


                guest_atom_pos_list = [atom.pos for atom in guest_compound.atoms()]

                atom_indicies = np.where(np.logical_or(np.any(guest_atom_pos_list < host_bounds.mins, axis=1), np.any(guest_atom_pos_list > host_bounds.maxes, axis=1)))[0]
                for ai in atom_indicies:
                    # print("Guest atom {} outside host box {} -- removing guest atom".format(guest_atoms[ai],host_bounds))
                    atoms_to_remove.add(guest_atoms[ai])
                    atoms_to_remove.update(guest_atoms[ai].bonded_atoms())

                guest.remove(atoms_to_remove)

                # remove overlapping atoms, and anything bonded to them (recursively)
                atoms_to_remove = set()
                for guest_atom in guest.atoms():
                    _, neighbors = kdtree.query(guest_atom.pos, k=10)
                    for host_atom_idx in neighbors:
                        if host_atom_idx < len(host_atom_list):
                            host_atom = host_atom_list[host_atom_idx]
                            if host_compound.min_periodic_distance(host_atom, guest_atom) < (vdw_radius(host_atom) + vdw_radius(guest_atom)):
                                # print("Guest atom {} overlaps with host atom {} -- removing guest atom".format(guest_atom,host_atom))
                                atoms_to_remove.add(guest_atom)
                                atoms_to_remove.update(guest_atom.bonded_atoms())

                guest.remove(atoms_to_remove)
                host_compound.add(guest, "guest_{}_{}_{}".format(xi,yi,zi))

if __name__ == "__main__":
    print "hello"



