import numpy as np
from scipy.spatial import cKDTree


def rank_points_by_density(points, k=14, box_min=None, box_max=None, wall_cutoff=0.0):
    """Rank points from lowest local density to highest.

    Notes
    -----
    If box_min, box_max and wall_cutoff are given, then points within wall_cutoff
    to the box boundaries are ignored. This ensure the effective "density" resulting
    from boundaries is acounted for.

    This can be useful for stringing together mutliple random walks
    with and without volume constraints. The use case for this function would be to find an existing site
    from a random walk, or another path, that has relatively lower density than other sites in the path.
    A new random walk path can begin from this site.

    If you want to find an unoccupied point with the lowest local density
    use mbuild.utils.density.find_low_density_point instead.

    Parameters
    ----------
    points : ndarray, shape (N, d)
        Coordinates of points.
    k : int, optional, default 14
        Number of neighbors to use for local density.
    box_min : array-like, shape (d,)
        Minimum boundary of box. Points closer than `wall_cutoff` will be ignored.
    box_max : array-like, shape (d,)
        Maximum boundary of box.
    wall_cutoff : float (nm)
        Distance from walls to ignore points.

    Returns
    -------
    sorted_indices : ndarray
        Indices of points (original array) sorted from lowest to highest local density.
    """
    points = np.asarray(points)
    N, dim = points.shape
    # Check for points too close to boundaries (if given)
    # Ignore these from density calculation, enforce high density due to proximity to wall
    if box_min is not None and box_max is not None and wall_cutoff > 0.0:
        box_min = np.asarray(box_min)
        box_max = np.asarray(box_max)
        mask = np.all(
            (points > box_min + wall_cutoff) & (points < box_max - wall_cutoff), axis=1
        )
        valid_indices = np.nonzero(mask)[0]
        filtered_points = points[mask]
    else:
        filtered_points = points
        valid_indices = np.arange(N)

    tree = cKDTree(filtered_points)
    dists, idxs = tree.query(filtered_points, k=k + 1)
    avg_dist = np.mean(dists[:, 1:], axis=1)
    local_density = 1 / (avg_dist + 1e-12)
    # Sort indices by increasing density
    sorted_order = np.argsort(local_density)
    return valid_indices[sorted_order]


def find_low_density_point(points, box_min, box_max, edge_buffer=0, n_candidates=5000):
    """Find an unoccupied point inside a box that is furthest from existing points.

    Parameters
    ----------
    points : ndarray, shape (N, 3)
        Array of existing points.
    box_min : array-like, shape (d, 3)
        Minimum coordinates of the box.
    box_max : array-like, shape (d, 3)
        Maximum coordinates of the box.
    edge_buffer : float (nm) default 0.0
        The buffer to prevent selection of coordinates within
        the buffer distance to the edge of the box.
    n_candidates : int
        Number of random candidate points to try.
        These points will be used to chose the point with lowest local density.
        Higher values may suffer performance costs, but give improved
        sampling.

    Returns
    -------
    best_point : ndarray, shape (d,3)
        Coordinates of the lowest-density point.
    """
    points = np.asarray(points)
    dim = points.shape[1]
    tree = cKDTree(points)
    # Create random candidates inside the box to test and sample from
    candidates = np.random.uniform(
        box_min + edge_buffer, box_max - edge_buffer, size=(n_candidates, dim)
    )
    dists, _ = tree.query(candidates, k=1)
    idx = np.argmax(dists)
    return candidates[idx]
