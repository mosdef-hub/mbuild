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


def apply_mask(host, guest, mask, guest_port_name="port"):
    """ """
    box = host.boundingbox(excludeG=False)

    mask = mask * box.lengths + box.mins

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
        equivalence_transform(brush, brush.labels[guest_port_name], closest_port)
        host.add(brush)
        port_pos[closest_point_idx,:] = np.array([np.inf, np.inf, np.inf])


def random_mask_3d(num_sites):
    """ """
    mask = np.random.random((num_sites, 3))
    return mask


def random_mask_2d(num_sites):
    """ """
    mask = random_mask_3d(num_sites)
    mask[:, 2] = 0
    return mask


def grid_mask_2d(n, m):
    """ """
    mask = np.zeros(shape=(n*m, 3), dtype=float)
    for i in range(n):
        for j in range(m):
            mask[i*m + j, 0] = i / n
            mask[i*m + j, 1] = j / m
    return mask


def grid_mask_3d(n, m, l):
    """ """
    mask = np.zeros(shape=(n*m*l, 3), dtype=float)
    for i in range(n):
        for j in range(m):
            for k in range(l):
                mask[i*m*l + j*l + k, 0] = i / n
                mask[i*m*l + j*l + k, 1] = j / m
                mask[i*m*l + j*l + k, 2] = k / l
    return mask


def sphere_mask(N):
    """Generate N evenly distributed points on the unit sphere.
    
    Sphere is centered at the origin. Alrgorithm based on the 'Golden Spiral'.

    Code by Chris Colbert from the numpy-discussion list:
    http://mail.scipy.org/pipermail/numpy-discussion/2009-July/043811.html

    """
    phi = (1 + np.sqrt(5)) / 2  # the golden ratio
    long_incr = 2*np.pi / phi   # how much to increment the longitude

    dz = 2.0 / float(N)         # a unit sphere has diameter 2
    bands = np.arange(N)        # each band will have one point placed on it
    z = bands * dz - 1 + (dz/2) # the height z of each band/point
    r = np.sqrt(1 - z*z)        # project onto xy-plane
    az = bands * long_incr      # azimuthal angle of point modulo 2 pi
    x = r * np.cos(az)
    y = r * np.sin(az)
    return np.column_stack((x, y, z))


def disk_mask(n):
    """ """
    radius = np.sqrt(np.arange(n) / float(n))

    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(n)

    points = np.zeros((n, 2))
    points[:,0] = np.cos(theta)
    points[:,1] = np.sin(theta)
    points *= radius.reshape((n, 1))

    return points


def vdw_radius(atomic_number):
    return 1.5


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

def add_bond(compound, type_A, type_B, dmin, dmax):
    """Ai-Bj distance is in [dmin, dmax] => add bond A1xB(Ai,Bj) (symmetric)."""
    for a1 in compound.atom_list_by_kind(type_A):
        nearest = compound.atoms_in_range(a1.pos, dmax, kind=type_B)
        for b1 in nearest:
            if (b1.kind==type_B) and (dmin <= compound.min_periodic_distance(b1.pos, a1.pos) <= dmax):
                compound.add(Bond(a1, b1))

if __name__ == "__main__":
    print "hello"



