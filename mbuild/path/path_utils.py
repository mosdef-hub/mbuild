"""Utility functions (mostly numba) for mbuild path generation.
These methods are primarly used by other functions and classes in mbuild.Path and
should rarely need to be called directly by users.
"""

import logging

import numpy as np
from numba import njit

logger = logging.getLogger(__name__)


@njit(cache=True, fastmath=True)
def random_coordinate(
    pos1,
    pos2,
    bond_length,
    thetas,
    r_vectors,
):
    """Default next_step method for HardSphereRandomWalk.
    This method takes in a a batch of thetas and vectors
    and creates a batch of random coordinates.

    Parameters
    ----------
    pos1 : np.ndarray (1,3), required
        The coordinate of the last accepted site.
    pos2 : np.ndarray (1,3), required
        The coordinate of the second to last accepted site.
    bond_length : float, required
        The fixed bond length between pos1 and the new coordinates.
    thetas : array-like (N, 1), required
        A set of possible angles used to determine new
        coordinates relative to pos2-pos1-new angle
    r_vectors : array-like (N, 3), required
        A set of normal vectors used perform rotations around.
    """

    if pos1 is None:  # pick random point in sphere.
        r_norm = compute_norms(r_vectors)
        return (pos2 + r_vectors / r_norm * bond_length).astype(np.float32)
    # pos1 and pos2 are defined, use available angles to sample new coordinates
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
def check_path(
    existing_points,
    new_point,
    radius,
    tolerance,
    pbc=np.array([False, False, False], dtype=bool),
    box_lengths=np.array([np.inf, np.inf, np.inf], dtype=np.float32),
    excluded_indices=np.empty(0, dtype=np.int64),
):
    """Default check path method for HardSphereRandomWalk.

    Parameters
    ----------
    existing_points : np.ndarray (N, 3), required
        Array of all fixed points in the system.
        This includes accepted sites from previous random walk steps
        and coordinates from included compounds.
    new_point : np.ndarray (1,3), required
        A candidate point for the next step.
    radius : float, required
        The radius used for hard-sphere overlap checks.
    tolerance : float, required
        Tolerance in center-to-center distances, allowing for rounding errors.
    pbc : np.ndarray (3,) of bool, optional
        Periodic boundary flags per axis. When an axis is periodic, the
        center-to-center distance along that axis uses the minimum-image
        convention so overlaps across a periodic face are detected. The
        default of all ``False`` reproduces the non-periodic behavior.
    box_lengths : np.ndarray (3,) of float, optional
        Box length along each axis, used for the minimum-image wrap on
        periodic axes. Only consulted where ``pbc`` is ``True``.
    excluded_indices : np.ndarray (K,) of int, optional
        Indices into ``existing_points`` that are skipped, sorted ascending.
        Used to leave out sites bonded to ``new_point``, whose separation is
        set by the bond length rather than by ``radius``. Indices at or past
        the end of ``existing_points`` are ignored. The default of an empty
        array checks every point.
    """
    if existing_points is None or existing_points.size == 0:
        return True
    min_sq_dist = (radius - tolerance) ** 2
    n_excluded = excluded_indices.shape[0]
    next_excluded = 0 # counter to increment along excluded_indices, which is sorted
    use_pbc = np.array([pbc[j] and not (box_lengths[j] + 1 == box_lengths[j]) for j in range(3)])
    for i in range(existing_points.shape[0]):
        # excluded_indices ascends alongside i, so one comparison per point
        if next_excluded < n_excluded and i == excluded_indices[next_excluded]:
            next_excluded += 1
            continue
        dist_sq = 0.0
        for j in range(existing_points.shape[1]):
            diff = existing_points[i, j] - new_point[j]
            # Apply minimum-image convention on periodic axes
            if use_pbc[j]:
                diff -= np.round(diff / box_lengths[j]) * box_lengths[j]
            dist_sq += diff * diff
        if dist_sq < min_sq_dist:
            return False
    return True


@njit(cache=True, fastmath=True)
def target_sq_distances(
    target_coordinate,
    new_points,
    pbc=np.array([False, False, False], dtype=bool),
    box_lengths=np.array([np.inf, np.inf, np.inf], dtype=np.float32),
):
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
def local_density(candidate, target_coords, r_cut):
    """Return number of target-type sites within r_cut of candidate."""
    r2_cut = r_cut * r_cut
    density = 0
    for i in range(target_coords.shape[0]):
        dx = candidate[0] - target_coords[i, 0]
        dy = candidate[1] - target_coords[i, 1]
        dz = candidate[2] - target_coords[i, 2]
        dist2 = dx * dx + dy * dy + dz * dz
        if dist2 < r2_cut:
            density += 1
    return density


@njit(cache=True, fastmath=True)
def target_density(candidates, target_coords, r_cut):
    """For a batch of candidate sites, calculate local density of target site-types."""
    n = candidates.shape[0]
    out = np.empty(n, dtype=np.float32)
    for i in range(n):
        out[i] = local_density(candidates[i], target_coords, r_cut)
    return out


@njit(cache=True, fastmath=True)
def norm(vec):
    """Use in place of np.linalg.norm inside of numba functions."""
    s = 0.0
    for i in range(vec.shape[0]):
        s += vec[i] * vec[i]
    return np.sqrt(s)


@njit
def compute_norms(vec):
    """Compute norms for multiple vectors."""
    n = vec.shape[0]
    r_norm = np.zeros((n, 1))
    for i in range(n):
        r_norm[i, 0] = norm(vec[i])  # Call your norm function
    return r_norm


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


@njit(cache=True, fastmath=True)
def calculate_sq_distances(
    target_coordinate,
    new_points,
    pbc=np.array([False, False, False], dtype=bool),
    box_lengths=np.array([np.inf, np.inf, np.inf], dtype=np.float32),
):
    """
    Return squared distances from target_coordinate to new_points.

    Parameters
    ----------
    target_coordinate : np.ndarray
        Single coordinate [x, y, z]
    new_points : np.ndarray
        Array of coordinates with shape (n, 3)
    pbc : np.ndarray
        Boolean array indicating periodic boundary conditions for each axis
    box_lengths : np.ndarray
        Box lengths for periodic boundary conditions

    Returns
    -------
    np.ndarray
        Squared distances from target to each point in new_points
    """
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
def find_candidates_within_radius(
    target_coordinate,
    candidate_points,
    radius,
    pbc=np.array([False, False, False], dtype=bool),
    box_lengths=np.array([np.inf, np.inf, np.inf], dtype=np.float32),
):
    """
    Find indices of candidate points within radius of target coordinate.

    Parameters
    ----------
    target_coordinate : np.ndarray
        Single coordinate [x, y, z]
    candidate_points : np.ndarray
        Array of coordinates with shape (n, 3)
    radius : float
        Search radius
    pbc : np.ndarray
        Boolean array indicating periodic boundary conditions for each axis
    box_lengths : np.ndarray
        Box lengths for periodic boundary conditions

    Returns
    -------
    np.ndarray
        Boolean mask of points within radius
    """
    n_points = candidate_points.shape[0]
    r2_cut = radius * radius
    within_radius = np.empty(n_points, dtype=np.bool_)

    for i in range(n_points):
        dx = target_coordinate[0] - candidate_points[i, 0]
        dy = target_coordinate[1] - candidate_points[i, 1]
        dz = target_coordinate[2] - candidate_points[i, 2]
        # Apply PBC per-axis
        if pbc[0]:
            dx -= np.round(dx / box_lengths[0]) * box_lengths[0]
        if pbc[1]:
            dy -= np.round(dy / box_lengths[1]) * box_lengths[1]
        if pbc[2]:
            dz -= np.round(dz / box_lengths[2]) * box_lengths[2]

        dist2 = dx * dx + dy * dy + dz * dz
        within_radius[i] = dist2 <= r2_cut

    return within_radius


def _angle_range(angles_sampler):
    """Return the smallest and largest angle a sampler can produce.

    Returns None for distributions with unbounded support.
    """
    if angles_sampler.distribution == "uniform":
        return (
            float(angles_sampler.kwargs["low"]),
            float(angles_sampler.kwargs["high"]),
        )
    if angles_sampler.distribution == "choice":
        angles = np.asarray(angles_sampler.kwargs["a"], dtype=float)
        return float(angles.min()), float(angles.max())
    return None


def check_angle_range(bond_length, radius, angles_sampler):
    """Compare the sampled angle range against the angle radius allows.

    Sites two bonds apart are separated by ``2 * bond_length * sin(theta / 2)``
    for a bond angle theta, and are not excluded from the overlap check. Angles
    below the critical angle place a new site within ``radius`` of the site two
    bonds back, so those angles are always rejected.

    Raises
    ------
    ValueError
        If no angle in the sampled range can avoid the overlap.
    """
    ratio = radius / (2.0 * bond_length)
    if ratio > 1.0:
        raise ValueError(
            f"A {radius=} larger than twice {bond_length=} leaves no bond angle "
            "that avoids overlapping the site two bonds back. Reduce radius or "
            "increase bond_length."
        )
    critical_angle = 2.0 * np.arcsin(ratio)
    angle_range = _angle_range(angles_sampler)
    if angle_range is None:
        return
    low, high = angle_range
    if critical_angle >= high:
        raise ValueError(
            f"With {bond_length=} and {radius=}, bond angles below "
            f"{np.degrees(critical_angle):.1f} degrees overlap the site two "
            f"bonds back, which rejects every angle in the sampled range of "
            f"{np.degrees(low):.1f} to {np.degrees(high):.1f} degrees. "
            "Reduce radius, increase bond_length, or raise rw_angles."
        )
    if critical_angle > low:
        if angles_sampler.distribution == "uniform":
            fraction = (critical_angle - low) / (high - low)
        else:
            angles = np.asarray(angles_sampler.kwargs["a"], dtype=float)
            fraction = float(np.mean(angles < critical_angle))
        logger.warning(
            f"With {bond_length=} and {radius=}, bond angles below "
            f"{np.degrees(critical_angle):.1f} degrees overlap the site two "
            f"bonds back. This rejects {fraction:.0%} of the sampled range "
            f"starting at {np.degrees(low):.1f} degrees, biasing the walk "
            "toward wider angles."
        )
