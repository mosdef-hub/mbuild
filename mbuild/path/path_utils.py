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
    return next_positions


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
def batch_rotate_molecule(
    coords, hinge_point, torsion_axis, bend_axes, bend_thetas, torsion_phis
):
    """Default next_step method for CompoundRandomWalk.

    Given a set of original molecular coordinates and a batch of bending axes, angles, and torsion angles
    a batch of updated trial molecular coordinates is created.
    """
    coords = coords.astype(np.float32)
    N_atoms = coords.shape[0]
    N_trials = bend_axes.shape[0]
    rotated_batch = np.zeros((N_trials, N_atoms, 3), dtype=np.float32)
    torsion_axis = torsion_axis / norm(torsion_axis)

    for t in range(N_trials):
        trial_coords = np.zeros((N_atoms, 3), dtype=np.float32)
        for i in range(N_atoms):
            trial_coords[i, :] = coords[i, :]
        for i in range(N_atoms):
            trial_coords[i, :] -= hinge_point
        # Apply torsion rotation before any bending.
        # Torsion done first as the torsion axis can change after bending.
        phi = torsion_phis[t]
        for i in range(N_atoms):
            trial_coords[i, :] = rotate_vector(trial_coords[i, :], torsion_axis, phi)
        # Apply bend rotation, bend_axis already normalized
        bend_axis = bend_axes[t] / norm(bend_axes[t])
        theta = bend_thetas[t]
        for i in range(N_atoms):
            # TODO: Do normalization here instead of batch generation?
            trial_coords[i, :] = rotate_vector(trial_coords[i, :], bend_axis, theta)
        # Reverse the original hinge point translation
        for i in range(N_atoms):
            trial_coords[i, :] += hinge_point
        rotated_batch[t, :, :] = trial_coords
    return rotated_batch


@njit(cache=True, fastmath=True)
def check_new_molecule(system_coords, new_coords, distance_tolerance):
    """Check if new molecule overlaps with system or with itself."""
    tol2 = distance_tolerance * distance_tolerance
    system_coords = system_coords.astype(np.float32)
    new_coords = new_coords.astype(np.float32)
    # Check for overlaps
    for i in range(new_coords.shape[0]):
        for j in range(system_coords.shape[0]):
            dx = new_coords[i, 0] - system_coords[j, 0]
            dy = new_coords[i, 1] - system_coords[j, 1]
            dz = new_coords[i, 2] - system_coords[j, 2]
            if dx * dx + dy * dy + dz * dz < tol2:
                return False
    return True


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
