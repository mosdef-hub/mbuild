from __future__ import division
import operator

from itertools import product

import numpy as np

from mbuild.coordinate_transform import (equivalence_transform, translate,
                                         rotate_around_x, rotate_around_y,
                                         rotate_around_z)
from mbuild.utils.validation import assert_port_exists
from mbuild import clone

__all__ = ['Pattern', 'DiskPattern', 'SpherePattern', 'Random2DPattern',
           'Random3DPattern', 'Grid2DPattern', 'Grid3DPattern']


class Pattern(object):
    def __init__(self, points, orientations=None):
        self.points = points
        if orientations is None:
            orientations = dict()
        self.orientations = orientations  # TODO: implement

    def __len__(self):
        return len(self.points)

    def __getitem__(self, item):
        return self.points[item]

    def scale(self, scalar):
        self.points *= scalar
        self._adjust_ports()

    def _adjust_ports(self):
        for orientation, ports in self.orientations.items():
            for port, point in zip(ports, self.points):
                translate(port, point)

    def apply(self, compound, orientation='', compound_port=''):
        """Arrange copies of a Compound as specified by the Pattern.

        Parameters
        ----------
        compound
        orientation

        Returns
        -------

        """
        compounds = list()
        if self.orientations.get(orientation):
            for port in self.orientations[orientation]:
                new_compound = clone(compound)
                new_port = new_compound.labels[compound_port]
                equivalence_transform(new_compound, new_port['up'], port['up'])

                compounds.append(new_compound)
        else:
            for point in self.points:
                new_compound = clone(compound)
                translate(new_compound, point)

                compounds.append(new_compound)
        return compounds

    def apply_to_compound(self, guest, guest_port_name='down', host=None,
                          backfill=None, backfill_port_name='up'):
        """Attach copies of a guest Compound to Ports on a host Compound.

        Parameters
        ----------
        guest
        guest_port_name
        host
        backfill
        backfill_port_name

        Returns
        -------

        """
        n_ports = len(host.available_ports())
        assert n_ports >= self.points.shape[0], "Not enough ports for pattern."

        assert_port_exists(guest_port_name, guest)
        box = host.boundingbox
        pattern = self.points * box.lengths + box.mins

        port_positions = np.empty(shape=(n_ports, 3))
        port_list = list()
        for port_idx, port in enumerate(host.available_ports()):
            port_positions[port_idx, :] = port['up']['middle'].pos
            port_list.append(port)

        used_ports = set()  # Keep track of used ports for backfilling.
        guests = []
        for point in pattern:
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


class Random2DPattern(Pattern):
    def __init__(self, n, orientations=None):
        points = np.random.random((n, 3))
        points[:, 2] = 0
        super(Random2DPattern, self).__init__(points=points, orientations=orientations)


class Random3DPattern(Pattern):
    def __init__(self, n, orientations=None):
        points = np.random.random((n, 3))
        super(Random3DPattern, self).__init__(points=points, orientations=orientations)


class Grid2DPattern(Pattern):
    def __init__(self, n, m, orientations=None):
        points = np.zeros(shape=(n*m, 3), dtype=float)
        for i, j in product(range(n), range(m)):
            points[i*m + j, 0] = i / n
            points[i*m + j, 1] = j / m
        super(Grid2DPattern, self).__init__(points=points, orientations=orientations)


class Grid3DPattern(Pattern):
    def __init__(self, n, m, l, orientations=None):
        points = np.zeros(shape=(n*m*l, 3), dtype=float)
        for i, j, k in product(range(n), range(m), range(l)):
            points[i*m*l + j*l + k, 0] = i / n
            points[i*m*l + j*l + k, 1] = j / m
            points[i*m*l + j*l + k, 2] = k / l
        super(Grid3DPattern, self).__init__(points=points, orientations=orientations)


class SpherePattern(Pattern):
    """Generate N evenly distributed points on the unit sphere.

    Sphere is centered at the origin. Alrgorithm based on the 'Golden Spiral'.

    Code by Chris Colbert from the numpy-discussion list:
    http://mail.scipy.org/pipermail/numpy-discussion/2009-July/043811.html

    """
    def __init__(self, n):
        phi = (1 + np.sqrt(5)) / 2  # the golden ratio
        long_incr = 2*np.pi / phi   # how much to increment the longitude

        dz = 2.0 / float(n)         # a unit sphere has diameter 2
        bands = np.arange(n)        # each band will have one point placed on it
        z = bands * dz - 1 + (dz/2) # the height z of each band/point
        r = np.sqrt(1 - z*z)        # project onto xy-plane
        az = bands * long_incr      # azimuthal angle of point modulo 2 pi
        x = r * np.cos(az)
        y = r * np.sin(az)
        points = np.column_stack((x, y, z))

        from mbuild.port import Port
        ports = list()
        for point in points:
            port = Port()
            ports.append(port)
            # Make the top of the port point toward the positive x axis.
            rotate_around_z(port, -np.pi/2)
            # Raise up (or down) the top of the port in the z direction.
            rotate_around_y(port, -np.arcsin(point[2]))
            # Rotate the Port along the z axis.
            rotate_around_z(port, np.arctan2(point[1], point[0]))
            # Move the Port a bit away from the surface of the Sphere.
            #translate(port, point + 0.07)

        super(SpherePattern, self).__init__(points=points,
                                            orientations={'normal': ports})


class DiskPattern(Pattern):
    """ """
    def __init__(self, n, orientations=None):
        radius = np.sqrt(np.arange(n) / float(n))

        golden_angle = np.pi * (3 - np.sqrt(5))
        theta = golden_angle * np.arange(n)

        points = np.zeros((n, 3))
        points[:, 0] = np.cos(theta)
        points[:, 1] = np.sin(theta)
        points *= radius.reshape((n, 1))
        super(DiskPattern, self).__init__(points=points, orientations=orientations)
