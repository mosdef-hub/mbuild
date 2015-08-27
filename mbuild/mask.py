from __future__ import division

from itertools import product

import numpy as np

from mbuild.coordinate_transform import equivalence_transform
from mbuild.utils.validation import assert_port_exists
from mbuild import clone

__all__ = ['apply_mask', 'random_mask_2d', 'random_mask_3d', 'sphere_mask',
           'grid_mask_2d', 'grid_mask_3d', 'disk_mask']

def apply_mask(host, guest, mask, guest_port_name='down', backfill=None,
               backfill_port_name='up'):
    """Attach guest Compounds to a host Compound in the pattern of a mask.

    Args:
        host (mbuild.Compound):
        guest (mbuild.Compound):
        mask (np.ndarray):
        guest_port_name (str):
        backfill (Compound, optional):
    """
    assert_port_exists(guest_port_name, guest)
    box = host.boundingbox
    mask = mask * box.lengths + box.mins

    n_ports = len(host.referenced_ports())
    assert n_ports >= mask.shape[0], "Not enough ports for mask."

    port_positions = np.empty(shape=(n_ports, 3))
    port_list = list()
    for port_idx, port in enumerate(host.referenced_ports()):
        port_positions[port_idx, :] = port.up.middle.pos
        port_list.append(port)

    used_ports = set()  # Keep track of used ports for backfilling.
    guests = []
    for point in mask:
        closest_point_idx = np.argmin(host.min_periodic_distance(point, port_positions))
        closest_port = port_list[closest_point_idx]
        used_ports.add(closest_port)

        # Attach the guest to the closest port.
        new_guest = clone(guest)
        equivalence_transform(new_guest, new_guest.labels[guest_port_name], closest_port)
        guests.append(new_guest)

        # Move the port as far away as possible (simpler than removing it).
        # There may well be a more elegant/efficient way of doing this.
        port_positions[closest_point_idx, :] = np.array([np.inf, np.inf, np.inf])

    backfills = []
    if backfill:
        assert_port_exists(backfill_port_name, backfill)
        # Attach the backfilling Compound to unused ports.
        for port in port_list:
            if port not in used_ports:
                new_backfill = clone(backfill)
                # Might make sense to have a backfill_port_name option...
                equivalence_transform(
                    new_backfill, new_backfill.labels[backfill_port_name], port)
                backfills.append(new_backfill)
    return guests, backfills


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
    for i, j in product(range(n), range(m)):
        mask[i*m + j, 0] = i / n
        mask[i*m + j, 1] = j / m
    return mask


def grid_mask_3d(n, m, l):
    """ """
    mask = np.zeros(shape=(n*m*l, 3), dtype=float)
    for i, j, k in product(range(n), range(m), range(l)):
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
    points[:, 0] = np.cos(theta)
    points[:, 1] = np.sin(theta)
    points *= radius.reshape((n, 1))

    return points
