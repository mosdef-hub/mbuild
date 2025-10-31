import numpy as np
from numba import njit


class Constraint:
    def __init__(self):
        pass

    def is_inside(self, points, buffer):
        """Return True if point satisfies constraint (inside), else False."""
        raise NotImplementedError("Must be implemented in subclasses")


class CuboidConstraint(Constraint):
    def __init__(self, Lx, Ly, Lz, center=(0, 0, 0)):
        self.center = np.asarray(center)
        self.mins = self.center - np.array([Lx / 2, Ly / 2, Lz / 2])
        self.maxs = self.center + np.array([Lx / 2, Ly / 2, Lz / 2])

    def is_inside(self, points, buffer):
        """Points and buffer are passed in from HardSphereRandomWalk"""
        return is_inside_cuboid(
            mins=self.mins, maxs=self.maxs, points=points, buffer=buffer
        )


class SphereConstraint(Constraint):
    def __init__(self, center, radius):
        self.center = np.array(center)
        self.radius = radius
        self.mins = self.center - self.radius
        self.maxs = self.center + self.radius

    def is_inside(self, points, buffer):
        """Points and buffer are passed in from HardSphereRandomWalk"""
        return is_inside_sphere(points=points, sphere_radius=self.radius, buffer=buffer)


class CylinderConstraint(Constraint):
    def __init__(self, center, radius, height):
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
        """Points and buffer are passed in from HardSphereRandomWalk"""
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
