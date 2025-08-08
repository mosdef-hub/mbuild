import numpy as np
from numba import njit


class Constraint:
    def __init__(self):
        pass

    def is_inside(self, point):
        """Return True if point satisfies constraint (inside), else False."""
        raise NotImplementedError("Must be implemented in subclasses")


class CuboidConstraint(Constraint):
    def __init__(self, Lx, Ly, Lz, center=(0, 0, 0)):
        self.center = np.asarray(center)
        self.mins = self.center - np.array([Lx, Ly, Lz])
        self.maxs = self.center + np.array([Lx, Ly, Lz])

    # Point to actual method to use
    def is_inside(self, points, radius):
        return is_inside_cuboid(
            mins=self.mins, maxs=self.maxs, points=points, radius=radius
        )


class SphereConstraint(Constraint):
    def __init__(self, center, radius):
        self.center = np.array(center)
        self.radius = radius

    def is_inside(self, points, radius):
        return is_inside_sphere(points=points, sphere_radius=self.radius, radius=radius)


@njit(cache=True, fastmath=True)
def is_inside_sphere(sphere_radius, points, radius):
    n_points = points.shape[0]
    results = np.empty(n_points, dtype=np.bool_)
    max_distance = sphere_radius - radius
    max_distance_sq = max_distance * max_distance
    for i in range(n_points):
        dist_from_center_sq = 0.0
        for j in range(3):
            dist_from_center_sq += points[i, j] * points[i, j]
        results[i] = dist_from_center_sq < max_distance_sq
    return results


@njit(cache=True, fastmath=True)
def is_inside_cuboid(mins, maxs, points, radius):
    n_points = points.shape[0]
    results = np.empty(n_points, dtype=np.bool_)
    for i in range(n_points):
        inside = True
        for j in range(3):
            if points[i, j] - radius < mins[j] or points[i, j] + radius > maxs[j]:
                inside = False
                break
        results[i] = inside
    return results
