"""Utility functions (mostly numba) for mbuild path generation"""

import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def random_coordinate(
    pos1,
    pos2,
    bond_length,
    thetas,
    r_vectors,
    batch_size,
):
    """Default next_step method for HardSphereRandomWalk."""
    v1 = pos2 - pos1
    v1_norm = v1 / norm(v1)
    dot_products = (r_vectors * v1_norm).sum(axis=1)
    r_perp = r_vectors - dot_products[:, None] * v1_norm
    norms = np.sqrt((r_perp * r_perp).sum(axis=1))
    # Handle rare cases where rprep vectors approach zero
    norms = np.where(norms < 1e-6, 1.0, norms)
    r_perp_norm = r_perp / norms[:, None]
    # Batch of trial next-step vectors using angles and r_norms
    cos_thetas = np.cos(thetas)
    sin_thetas = np.sin(thetas)
    v2s = cos_thetas[:, None] * v1_norm + sin_thetas[:, None] * r_perp_norm
    # Batch of trial positions
    next_positions = pos1 + v2s * bond_length
    return next_positions.astype(np.float32)


@njit(cache=True, fastmath=True)
def check_path(existing_points, new_point, radius, tolerance):
    """Default check path method for HardSphereRandomWalk."""
    min_sq_dist = (radius - tolerance) ** 2
    for i in range(existing_points.shape[0]):
        dist_sq = 0.0
        for j in range(existing_points.shape[1]):
            diff = existing_points[i, j] - new_point[j]
            dist_sq += diff * diff
        if dist_sq < min_sq_dist:
            return False
    return True


@njit(cache=True, fastmath=True)
def target_sq_distances(target_coordinate, new_points, pbc=[False, False, False], box_lengths=[None, None, None]):
    """Return squared distances from target_coordinate to new_points."""
    n_points = new_points.shape[0]
    sq_distances = np.empty(n_points, dtype=np.float32)
    for i in range(n_points):
        dx = target_coordinate[0] - new_points[i, 0]
        dy = target_coordinate[1] - new_points[i, 1]
        dz = target_coordinate[2] - new_points[i, 2]
        # Apply PBC per-axis
        if pbc[0]:
            dx -= np.round(dx / box_lengths[0]) * box_lengths[0]
        if pbc[1]:
            dy -= np.round(dy / box_lengths[1]) * box_lengths[1]
        if pbc[2]:
            dz -= np.round(dz / box_lengths[2]) * box_lengths[2]
        sq_distances[i] = dx * dx + dy * dy + dz * dz
    return sq_distances


@njit(cache=True, fastmath=True)
def norm(vec):
    """Use in place of np.linalg.norm inside of numba functions."""
    s = 0.0
    for i in range(vec.shape[0]):
        s += vec[i] * vec[i]
    return np.sqrt(s)


@njit(cache=True, fastmath=True)
def rotate_vector(v, axis, theta):
    """Rotate vector v around a normalized axis by angle theta using Rodrigues' formula."""
    c = np.cos(theta)
    s = np.sin(theta)
    k = axis
    k_dot_v = k[0] * v[0] + k[1] * v[1] + k[2] * v[2]
    cross = np.zeros(3)
    cross[0] = k[1] * v[2] - k[2] * v[1]
    cross[1] = k[2] * v[0] - k[0] * v[2]
    cross[2] = k[0] * v[1] - k[1] * v[0]
    rotated = np.zeros(3)
    for i in range(3):
        rotated[i] = v[i] * c + cross[i] * s + k[i] * k_dot_v * (1 - c)
    return rotated
