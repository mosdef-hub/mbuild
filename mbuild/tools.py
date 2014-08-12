from __future__ import division
from copy import deepcopy
from atom import Atom
from compound import Compound
from mbuild.coordinate_transform import equivalence_transform

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


def solvate(host_compound, guest_compound, host_bounds, guest_bounds):
    assert(np.shape())
    # we may want to make sure that the axes of the two boxes line up

    # replicate the quest so that it's bigger than the host



    host_compound.add(guest_compound, 'guest')



    for ag in guest_compound.atoms():
        for neighbor in ag.getAtomsInRange(o.pos, 10, maxItems=10, kind='*'):
            if neighbor is o:
                continue
            neighbor_sigma = Prototype.getAttr(neighbor.kind, "sigma", default=float("-inf"))
            if isinstance(neighbor_sigma, Quantity):
                neighbor_sigma = neighbor_sigma._value
            print "neighbor_sigma=" + str(neighbor_sigma)
            if o.distance(neighbor) < o_sigma + neighbor_sigma:
                print str(o) + " is close to " + str(neighbor)
                print "removing " + str(o)

                for b in o.bonds:
                    print str(b)
                    self.remove(b.atom1)
                    self.remove(b.atom2)
                    print "removing " + str(b.atom1)
                    print "removing " + str(b.atom2)






if __name__ == "__main__":
    print "hello"



