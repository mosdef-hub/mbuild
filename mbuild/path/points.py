import logging

import numpy as np

from mbuild.exceptions import PathConvergenceError
from mbuild.path.constraints import CuboidConstraint, CylinderConstraint

logger = logging.getLogger(__name__)


def get_second_point(state, existing_points, beads, check_path, next_step):
    """Generate a secound point from the given first point using RandomWalkState.

    Candidates are generated around the chain tip using a batch of trial angles
    and vectors. If the walk uses ``link-linear`` connectivity and a previous
    direction is available, candidates are angle-constrained relative to that
    direction; otherwise a sphere of candidates is generated around the tip.

    Parameters
    ----------
    state : RandomWalkState
        The current state of the random walk, containing radius, bond_length,
        tolerance, connectivity, volume_constraint, bias, and include_compound.
    existing_points : np.ndarray, shape (N, 3)
        Live coordinates accepted so far, sliced to ``coordinates[:state.count]``.
        Must contain at least one point (the initial point).
    beads : np.ndarray of str, shape (N,)
        Bead names corresponding to ``existing_points``, sliced to the same
        live length.
    check_path : callable
        Overlap-check function with signature
        ``check_path(existing_points, new_point, radius, tolerance) -> bool``.
    next_step : callable
        Coordinate-generation function with signature
        ``next_step(pos1, pos2, bond_length, thetas, r_vectors) -> np.ndarray``.
        Pass ``pos1=None`` to generate a sphere around ``pos2``.

    Returns
    -------
    np.ndarray, shape (3,) or None
        The accepted second coordinate, or ``None`` if no valid candidate was
        found within the trial batch.

    """
    batch_angles, batch_vectors = generate_trials(state)
    # If this RW is using link linear, pos2 = last site of last
    # Set pos1 and pos2 before checking include compound and combining coordinates
    if state.connectivity == "link-linear" and len(existing_points) > 1:
        pos1 = existing_points[-1]
        pos2 = existing_points[-2]
    else:
        pos1 = None
        pos2 = existing_points[-1]
    # Update existing points to include those in the compound.
    if state.include_compound:
        existing_points = np.concat((existing_points, state.include_compound.xyz))
    xyzs = next_step(
        pos1=pos1,
        pos2=pos2,
        bond_length=state.bond_length,
        thetas=batch_angles,
        r_vectors=batch_vectors,
    )
    if state.volume_constraint:
        is_inside_mask = state.volume_constraint.is_inside(
            points=xyzs, buffer=state.radius
        )
        xyzs = xyzs[is_inside_mask]

    if state.bias:
        xyzs = state.bias(candidates=xyzs, coordinates=existing_points, names=beads)

    for xyz in xyzs:
        if check_path(
            existing_points=existing_points,
            new_point=xyz,
            radius=state.radius,
            tolerance=state.tolerance,
        ):
            return xyz
    return None


def get_initial_point(state, existing_points, beads, check_path, next_step):
    """Generate a starting pointt from a RandomWalkState.

    The strategy for choosing a starting point depends on ``state.initial_point``:

    - **np.ndarray (3,)**: the array is used directly as the starting coordinate.
    - **int**: treated as an index into ``existing_points``; a new point is
      generated in a sphere around that coordinate, filtered by volume
      constraint and bias, and checked for overlaps.
    - **None with volume_constraint**: candidates are sampled from the volume
      constraint's low-density regions and checked for overlaps.
    - **None without volume_constraint**: candidates are drawn uniformly at
      random within the bounding box of ``existing_points`` (or a unit sphere
      around the origin if no points exist yet) and checked for overlaps.

    If ``state.include_compound`` is set, its coordinates are appended to
    ``existing_points`` before overlap checks so that the new point does not
    clash with the included compound's atoms.

    Parameters
    ----------
    state : RandomWalkState
        The current state of the random walk, containing initial_point, radius,
        bond_length, tolerance, connectivity, volume_constraint, bias,
        include_compound, and the trial-generation parameters.
    existing_points : np.ndarray, shape (N, 3)
        Live coordinates accepted so far, sliced to ``coordinates[:state.count]``.
        May be empty at the start of the first walk.
    beads : np.ndarray of str, shape (N,)
        Bead names corresponding to ``existing_points``, sliced to the same
        live length. Passed to the bias for sequence-aware scoring.
    check_path : callable
        Overlap-check function with signature
        ``check_path(existing_points, new_point, radius, tolerance) -> bool``.
    next_step : callable
        Coordinate-generation function with signature
        ``next_step(pos1, pos2, bond_length, thetas, r_vectors) -> np.ndarray``.
        Used only when ``state.initial_point`` is an int.

    Returns
    -------
    np.ndarray, shape (3,)
        The accepted starting coordinate.

    Raises
    ------
    ValueError
        If ``state.initial_point`` is an int that is out of bounds for
        ``existing_points``.
    PathConvergenceError
        If no valid starting point can be found within the trial batch,
        regardless of which strategy is used.
    """
    if state.include_compound:
        existing_points = np.concat((existing_points, state.include_compound.xyz))

    # An initial point was manuallyl given in hard_sphere_random_walk, use that.
    if isinstance(state.initial_point, np.ndarray) and state.initial_point.shape == (
        3,
    ):
        return state.initial_point

    # Passing in an index to specify an initial point from already defined set of coordinates
    elif isinstance(state.initial_point, int):
        if state.initial_point >= len(existing_points):
            raise ValueError(
                f"You passed a starting index of {state.initial_point} "
                f"but there are only {len(existing_points)} existing points in the path."
            )
        # generate point off of current path coordinates
        starting_xyz = existing_points[state.initial_point]
        batch_angles, batch_vectors = generate_trials(state)
        # TODO: If building from a path with coordinates, can we try to get both pos1 and pos2?
        xyzs = next_step(
            pos1=None,  # will generate sphere of points around pos2
            pos2=starting_xyz,
            bond_length=state.bond_length,
            thetas=batch_angles,
            r_vectors=batch_vectors,
        )
        if state.volume_constraint:
            is_inside_mask = state.volume_constraint.is_inside(
                points=xyzs, buffer=state.radius
            )
            xyzs = xyzs[is_inside_mask]

        if state.bias:
            xyzs = state.bias(candidates=xyzs, coordinates=existing_points, names=beads)

        # Set up PBC info from volume constraints
        if isinstance(state.volume_constraint, CuboidConstraint):
            pbc = state.volume_constraint.pbc
            box_lengths = state.volume_constraint.box_lengths.astype(np.float32)
        elif isinstance(state.volume_constraint, CylinderConstraint):
            pbc = (False, False, state.volume_constraint.periodic_height)
            box_lengths = np.array(
                [
                    state.volume_constraint.radius * 2,
                    state.volume_constraint.radius * 2,
                    state.volume_constraint.height,
                ]
            ).astype(np.float32)
        else:
            pbc = (None, None, None)
            box_lengths = (None, None, None)

        for i in range(len(xyzs)):
            xyz = xyzs[i]
            if any(pbc):
                xyz = state.volume_constraint.mins + np.mod(
                    xyz - state.volume_constraint.mins, box_lengths
                )
            if check_path(  # check for overlaps
                existing_points=existing_points,
                new_point=xyz,
                radius=state.radius,
                tolerance=state.tolerance,
            ):
                return xyz
        raise PathConvergenceError(
            f"Unable to find a starting point at {starting_xyz} "
            "without overlapping particles. "
            "Check your `initial_point` argument."
        )
    # No initial point given, but there is a volume constraint to sample starting points from
    # TODO: Use find_low_density_point here instead?
    elif state.volume_constraint:
        xyzs = state.volume_constraint.sample_candidates(
            points=existing_points, n_candidates=300, buffer=state.radius + 0.1
        )
        for xyz in xyzs:
            if check_path(
                existing_points=existing_points,
                new_point=xyz,
                radius=state.radius,
                tolerance=state.tolerance,
            ):
                return xyz
        raise PathConvergenceError(
            "Unable to find a starting point without overlapping particles. "
            "The density of the volume constraint may be too high."
        )
    else:  # completely random point
        if len(existing_points) == 0:
            max_dist = state.radius
            min_dist = -1 * state.radius
        else:
            # TODO: Update seed based on initial_point
            max_dist = np.max(existing_points, axis=0)
            min_dist = np.min(existing_points, axis=0)
        xyzs = state.rng.uniform(low=min_dist, high=max_dist, size=(300, 3))
        for xyz in xyzs:
            if check_path(
                existing_points=existing_points,
                new_point=xyz,
                radius=state.radius,
                tolerance=state.tolerance,
            ):
                return xyz
        raise PathConvergenceError(
            "Unable to find a starting point without overlapping particles. "
            "The density of the volume constraint may be too high."
        )


class AnglesSampler:
    """
    TODO:Allow for passing a specific 2D weighted pre-defined sample, as opposed to a numpy distribution.
    This would use np.random.choice instead as the sample method.
    NOTES
    -----
    "uniform" distribution should use 'low' and 'high' as kwargs.
    "normal" distribution should use 'loc' and 'scale' as kwargs.
    """

    def __init__(self, distributionStr, kwargs, seed):
        # Create a generator object for high-quality random numbers [9]
        self.rng = np.random.default_rng(seed)
        if distributionStr.lower() == "uniform":
            self.sampler = self.rng.uniform
            assert "low" in kwargs
            assert "high" in kwargs
        elif distributionStr.lower() == "normal":
            self.sampler = self.rng.normal
            assert "loc" in kwargs
            assert "scale" in kwargs
        elif distributionStr == "choice":
            self.sampler = self.rng.choice
            assert "a" in kwargs  # p is not required
        else:
            raise NotImplementedError(
                f"Sample Distribution {distributionStr} not supported."
            )
        self.kwargs = kwargs

    def sample(self, size=None):
        return self.sampler(size=size, **self.kwargs)


def generate_trials(state):
    """Use normal or uniform sampling on angles, uniform sampling on radius."""
    thetas = state.angles.sample(size=state.trial_batch_size).astype(np.float32)
    r = state.rng.uniform(-0.5, 0.5, size=(state.trial_batch_size, 3)).astype(
        np.float32
    )
    return thetas, r


def generate_crosslink_init_points():
    pass
