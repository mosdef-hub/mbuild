from __future__ import division

from itertools import product

import numpy as np

from mbuild.coordinate_transform import force_overlap
from mbuild.utils.validation import assert_port_exists
from mbuild import clone

__all__ = ['Pattern', 'DiskPattern', 'SpherePattern', 'Random2DPattern',
           'Random3DPattern', 'Grid2DPattern', 'Grid3DPattern']


class Pattern(object):
    def __init__(self, points, orientations=None, scale=None, **kwargs):
        self.points = points
        if orientations is None:
            orientations = dict()
        self.orientations = orientations  # TODO: implement for more patterns
        if scale is not None:
            self.scale(scale)

    def __len__(self):
        return len(self.points)

    def __getitem__(self, item):
        return self.points[item]

    def scale(self, by):
        """Scale the points in the Pattern.

        Parameters
        ----------
        by : float or np.ndarray, shape=(3,)
            The factor to scale by. If a scalar, scale all directions isotropically.
            If np.ndarray, scale each direction independently.
        """
        self.points *= np.asarray([by])
        self._adjust_ports()

    def _adjust_ports(self):
        for orientation, ports in self.orientations.items():
            for port, point in zip(ports, self.points):
                port.translate(point)

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
                (new_compound, new_port['up'], port['up'])

                compounds.append(new_compound)
        else:
            for point in self.points:
                new_compound = clone(compound)
                new_compound.translate(point)

                compounds.append(new_compound)
        return compounds

    def apply_to_compound(self, guest, guest_port_name='down', host=None,
                          backfill=None, backfill_port_name='up', scale=True):
        """Attach copies of a guest Compound to Ports on a host Compound.

        Parameters
        ----------
        guest : mb.Compound
            The Compound prototype to be applied to the host Compound
        guest_port_name : str, optional, default='down'
            The name of the port located on `guest` to attach to the host
        host : mb.Compound, optional, default=None
            A Compound with available ports to add copies of `guest` to
        backfill : mb.Compound, optional, default=None
            A Compound to add to the remaining available ports on `host`
            after clones of `guest` have been added for each point in the
            pattern
        backfill_port_name : str, optional, default='up'
            The name of the port located on `backfill` to attach to the host
        scale : bool, optional, default=True
            Scale the points in the pattern to the lengths of the `host`'s
            `boundingbox` and shift them by the `boundingbox`'s mins

        Returns
        -------

        """
        n_ports = len(host.available_ports())
        assert n_ports >= self.points.shape[0], "Not enough ports for pattern."
        assert_port_exists(guest_port_name, guest)
        box = host.boundingbox
        if scale:
            self.scale(box.lengths)
            self.points += box.mins
        pattern = self.points
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
            force_overlap(new_guest, new_guest.labels[guest_port_name], closest_port)
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
                    force_overlap(new_backfill,
                                  new_backfill.labels[backfill_port_name],
                                  port)
                    backfills.append(new_backfill)
        return guests, backfills


class Random2DPattern(Pattern):
    def __init__(self, n, seed=None, **kwargs):
        """ Generate n random points on a 2D grid along z = 0

        Attributes
        ----------
        n : int
            Number of points to generate
        seed : int
            Seed for random number generation

        """

        if seed:
            np.random.seed(seed)
        points = np.random.random((n, 3))
        points[:, 2] = 0
        super(Random2DPattern, self).__init__(points=points, **kwargs)


class Random3DPattern(Pattern):
    def __init__(self, n, seed=None, **kwargs):
        """ Generate n random points on a 3D grid 

        Attributes
        ----------
        n : int
            Number of points to generate
        seed : int
            Seed for random number generation

        """

        if seed:
            np.random.seed(seed)
        points = np.random.random((n, 3))
        super(Random3DPattern, self).__init__(points=points, **kwargs)


class Grid2DPattern(Pattern):
    def __init__(self, n, m, **kwargs):
        """ Generate a 2D grid (n x m) of points along z = 0

        Notes
        -----
        Points span [0,1) along x and y axes

        Attributes
        ---------
        n : int
            Number of grid rows
        m : int
            Number of grid columns

        """

        points = np.zeros(shape=(n*m, 3), dtype=float)
        for i, j in product(range(n), range(m)):
            points[i*m + j, 0] = i / n
            points[i*m + j, 1] = j / m
        super(Grid2DPattern, self).__init__(points=points, **kwargs)


class Grid3DPattern(Pattern):
    def __init__(self, n, m, l, **kwargs):
        """ Generate a 3D grid (n x m x l) of points

        Notes
        -----
        Points span [0,1) along x, y, and z axes

        Attributes
        ---------
        n : int
            Number of grid rows
        m : int
            Number of grid columns
        l : int
            Number of grid aisles

        """

        points = np.zeros(shape=(n*m*l, 3), dtype=float)
        for i, j, k in product(range(n), range(m), range(l)):
            points[i*m*l + j*l + k, 0] = i / n
            points[i*m*l + j*l + k, 1] = j / m
            points[i*m*l + j*l + k, 2] = k / l
        super(Grid3DPattern, self).__init__(points=points, **kwargs)


class SpherePattern(Pattern):
    """Generate N evenly distributed points on the unit sphere.

    Sphere is centered at the origin. Alrgorithm based on the 'Golden Spiral'.

    Code by Chris Colbert from the numpy-discussion list:
    http://mail.scipy.org/pipermail/numpy-discussion/2009-July/043811.html

    """
    def __init__(self, n, **kwargs):
        phi = (1 + np.sqrt(5)) / 2   # the golden ratio
        long_incr = 2*np.pi / phi    # how much to increment the longitude
        dz = 2.0 / float(n)          # a unit sphere has diameter 2
        bands = np.arange(n)         # each band will have one point placed on it
        z = bands * dz - 1.0 + (dz/2.0)  # the height z of each band/point
        r = np.sqrt(1.0 - z*z)         # project onto xy-plane
        az = bands * long_incr       # azimuthal angle of point modulo 2 pi
        x = r * np.cos(az)
        y = r * np.sin(az)
        points = np.column_stack((x, y, z))

        from mbuild.port import Port
        if kwargs.get('orientations') is None:
            ports = list()
            for point in points:
                port = Port()
                ports.append(port)
                # Make the top of the port point toward the positive x axis.
                port.spin(-np.pi/2, [0, 0, 1])
                # Raise up (or down) the top of the port in the z direction.
                port.spin(-np.arcsin(point[2]), [0, 1, 0])
                # Rotate the Port along the z axis.
                port.spin(np.arctan2(point[1], point[0]), [0, 0, 1])
                # Move the Port a bit away from the surface of the Sphere.
                # translate(port, point + 0.07)
            kwargs['orientations'] = {'normal': ports}
        else:
            raise NotImplementedError('Custom orientation support is not yet '
                                      'implemented.')
        super(SpherePattern, self).__init__(points=points, **kwargs)


class DiskPattern(Pattern):
    """Generate N evenly distributed points on the unit circle along z = 0.

    Disk is centered at the origin. Algorithm based on Vogel's method.

    Code by Alexandre Devert:
    http://blog.marmakoide.org/?p=1

    """

    def __init__(self, n, **kwargs):
        radius = np.sqrt(np.arange(n) / float(n))
        golden_angle = np.pi * (3 - np.sqrt(5))
        theta = golden_angle * np.arange(n)
        points = np.zeros((n, 3))
        points[:, 0] = np.cos(theta)
        points[:, 1] = np.sin(theta)
        points *= radius.reshape((n, 1))
        super(DiskPattern, self).__init__(points=points, **kwargs)
