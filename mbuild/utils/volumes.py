import numpy as np
from numba import njit
from scipy.spatial import cKDTree


class Constraint:
    """
    Defines a volume that acts as a constraint in mbuild.path.HardSphereRandomWalk.
    This is the base class from which all constraints inherit from.

    Notes
    -----
    Design and implement your own volume constraint by inheriting from this class
    and implementing the constraint check in an `is_inside` method.

    `is_inside` is expected to return a mask of booleans.

    The existing contraints of CuboidConstraint, SphereConstraint, and CylinderConstraint
    call numba methods from their respectice `is_inside` method, but that is not a required
    implementation to design your own constraint.
    """

    def __init__(self):
        pass

    def is_inside(self, points, buffer):
        """Return True if point satisfies constraint (inside), else False."""
        raise NotImplementedError("Must be implemented in subclasses")

    def sample_candidates(self, points, n_candiates, buffer):
        """Sample the volume for canidate points sorted by lowest local density."""
        raise NotImplementedError("Must be implemented in subclasses")

    def find_low_density_points(self, points, n_candidates, buffer):
        low_density_points = self.sample_candidates(points=points, n_candidates=n_candidates, buffer=buffer)


class CuboidConstraint(Constraint):
    """Creates a cuboid constraint.

    Parameters
    ----------
    Lx : float, required
        The length of the volume along the x-axis.
    Ly : float, required
        The length of the volume along the y-axis.
    Lz : float, required
        The length of the volume along the z-axis.
    center : array-like (1,3), default = (0, 0, 0)
        Defines the center of the volume.
    """

    def __init__(self, Lx, Ly, Lz, center=(0, 0, 0)):
        self.center = np.asarray(center)
        self.mins = self.center - np.array([Lx / 2, Ly / 2, Lz / 2])
        self.maxs = self.center + np.array([Lx / 2, Ly / 2, Lz / 2])

    def is_inside(self, points, buffer):
        """Check a set of coordinates against the volume constraint.

        Parameters
        ----------
        points : ndarray (N, 3), required
            The set of points to check against the volume constraint.
        buffer : float, required
            Buffer used for rounding

        Returns
        -------
        Mask of booleans of length N corresponding to each point.
        """
        return is_inside_cuboid(
            mins=self.mins, maxs=self.maxs, points=points, buffer=buffer
        )

    def sample_candidates(self, points, n_candidates, buffer):
        """"""
        points = np.asarray(points)
        # No points to check against, randomly sample inside a cube
        if len(points) == 0:
            return np.random.uniform(self.mins + buffer, self.maxs - buffer, size=(n_candidates, 3))

        # Create random candidates inside the box to test and sample from
        candidates = np.random.uniform(
            self.mins + buffer, self.maxs - buffer, size=(n_candidates, 3)
        )
        tree = cKDTree(points)
        dists, _ = tree.query(candidates, k=1)
        sorted_order = np.argsort(-dists)
        return candidates[sorted_order]


class SphereConstraint(Constraint):
    """Creates a spherical constraint.

    Parameters
    ----------
    radius : float, required
        The radius of the sphere
    center : array-like (1,3), default = (0, 0, 0)
        Defines the center point of the sphere.
    """

    def __init__(self, center, radius):
        self.center = np.array(center)
        self.radius = radius
        self.mins = self.center - self.radius
        self.maxs = self.center + self.radius

    def is_inside(self, points, buffer):
        """Check a set of coordinates against the volume constraint.

        Parameters
        ----------
        points : ndarray (N, 3), required
            The set of points to check against the volume constraint.
        buffer : float, required
            Buffer used for rounding

        Returns
        -------
        Mask of booleans of length N corresponding to each point.
        """
        return is_inside_sphere(points=points, sphere_radius=self.radius, buffer=buffer)

    def sample_candidates(self, points, n_candidates, buffer):
        effective_radius = self.radius - buffer
        # Create random candidates inside the box to test and sample from
        dirs = np.random.normal(size=(n_candidates, 3))
        dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        u = np.random.random(size=n_candidates)
        radii = effective_radius * (u ** (1/3))
        # generate candidate points
        candidates = self.center + dirs * radii[:, None]
        # If just sampling from the volume, no KDTree needed
        if points is None or len(points) == 0:
            return candidates
        # Want to sample around existing points, sort by local density
        points = np.asarray(points)
        tree = cKDTree(points)
        dists, _ = tree.query(candidates, k=1)
        sorted_order = np.argsort(-dists)
        return candidates[sorted_order]


class CylinderConstraint(Constraint):
    """Creates a cylindrical constraint.

    Parameters
    ----------
    radius : float, required
        The radius of the cylinder
    height : float, required
        The height  of the cylinder
    center : array-like (1,3), default = (0, 0, 0)
        Defines the center point of the sphere.
    """

    def __init__(self, radius, height, center=(0, 0, 0)):
        self.center = np.array(center)
        self.height = height
        self.radius = radius
        self.mins = np.array(
            [
                self.center[0] - self.radius,
                self.center[1] - self.radius,
                self.center[2] - self.height / 2,
            ]
        )
        self.maxs = np.array(
            [
                self.center[0] + self.radius,
                self.center[1] + self.radius,
                self.center[2] + self.height / 2,
            ]
        )

    def is_inside(self, points, buffer):
        """Check a set of coordinates against the volume constraint.

        Parameters
        ----------
        points : ndarray (N, 3), required
            The set of points to check against the volume constraint.
        buffer : float, required
            Buffer used for rounding

        Returns
        -------
        Mask of booleans of length N corresponding to each point.
        """
        return is_inside_cylinder(
            points=points,
            center=self.center,
            cylinder_radius=self.radius,
            height=self.height,
            buffer=buffer,
        )


@njit(cache=True, fastmath=True)
def is_inside_cylinder(points, center, cylinder_radius, height, buffer):
    n_points = points.shape[0]
    results = np.empty(n_points, dtype=np.bool_)
    max_r = cylinder_radius - buffer
    max_r_sq = max_r * max_r
    half_height = height / 2.0
    max_z = half_height - buffer
    for i in range(n_points):
        dx = points[i, 0] - center[0]
        dy = points[i, 1] - center[1]
        dz = points[i, 2] - center[2]
        r_sq = dx * dx + dy * dy
        inside_radial = r_sq <= max_r_sq
        inside_z = abs(dz) <= max_z
        results[i] = inside_radial and inside_z
    return results


@njit(cache=True, fastmath=True)
def is_inside_sphere(sphere_radius, points, buffer):
    n_points = points.shape[0]
    results = np.empty(n_points, dtype=np.bool_)
    max_distance = sphere_radius - buffer
    max_distance_sq = max_distance * max_distance
    for i in range(n_points):
        dist_from_center_sq = 0.0
        for j in range(3):
            dist_from_center_sq += points[i, j] * points[i, j]
        results[i] = dist_from_center_sq < max_distance_sq
    return results


@njit(cache=True, fastmath=True)
def is_inside_cuboid(mins, maxs, points, buffer):
    n_points = points.shape[0]
    results = np.empty(n_points, dtype=np.bool_)
    for i in range(n_points):
        inside = True
        for j in range(3):
            if points[i, j] - buffer < mins[j] or points[i, j] + buffer > maxs[j]:
                inside = False
                break
        results[i] = inside
    return results
