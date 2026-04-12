import logging

import logging

import numpy as np

from mbuild.exceptions import PathConvergenceError
from mbuild.path.constraints import CuboidConstraint, CylinderConstraint

logger = logging.getLogger(__name__)



def get_second_point(state, existing_points, check_path, first_point):
    """Generate a secound point from the given first point using RandomWalkState."""
    # Find 2nd point
    found_valid_point = False
    if not state.volume_constraint:
        found_valid_point = False
        phiList = state.rng.uniform(0, 2 * np.pi, 10000)
        cos_thetaList = state.rng.uniform(-1, 1, 10000)
        thetaList = np.arccos(cos_thetaList)  # use to evenly sample sphereical space
        for phi, theta in zip(phiList, thetaList):
            offset = np.array(
                [
                    state.bond_length * np.sin(theta) * np.cos(phi),
                    state.bond_length * np.sin(theta) * np.sin(phi),
                    state.bond_length * np.cos(theta),
                ]
            )
            if check_path(  # check for overlaps
                existing_points=np.vstack((existing_points, first_point)),
                new_point=offset + first_point,
                radius=state.radius,
                tolerance=state.tolerance,
            ):
                found_valid_point = True
                break

    else:
        phiList = state.rng.uniform(0, 2 * np.pi, 10000)
        cos_thetaList = state.rng.uniform(-1, 1, 10000)
        thetaList = np.arccos(cos_thetaList)  # use to evenly sample sphereical space
        for phi, theta in zip(phiList, thetaList):
            offset = np.array(
                [
                    state.bond_length * np.sin(theta) * np.cos(phi),
                    state.bond_length * np.sin(theta) * np.sin(phi),
                    state.bond_length * np.cos(theta),
                ]
            )
            state.attempts += 1
            if state.termination.is_met() and not state.termination.success:
                logger.error("Random walk not successful.")
                logger.error(state.termination.summarize())
                return
            is_inside_mask = state.volume_constraint.is_inside(
                points=np.array([first_point + offset]), buffer=state.radius
            )
            if np.all(is_inside_mask) and check_path(  # check for overlaps
                existing_points=np.vstack((existing_points, first_point)),
                new_point=offset + first_point,
                radius=state.radius,
                tolerance=state.tolerance,
            ):
                found_valid_point = True
                break
    if not found_valid_point:
        raise PathConvergenceError(
            f"No viable second point found within constraint and next to {state.initial_point=}. Try using a smaller radius than {state.radius=}"
        )
    return first_point + offset


def get_initial_point(state, existing_points, check_path, next_step):
    """Generate a starting pointt from a RandomWalkState.
    Could initialize from:
        Specific xyz point
        Build off of an index in the current paths coordinates
        Random point in volume constraint
        Random point


    """
    if isinstance(state.initial_point, np.ndarray) and state.initial_point.shape == (
        3,
    ):
        return state.initial_point
    elif isinstance(state.initial_point, int) and state.initial_point < len(
        existing_points
    ):
        # generate point off of current path coordinates
        starting_xyz = existing_points[state.initial_point]
        batch_angles, batch_vectors = generate_trials(state)
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
            is_inside_mask = state.volume_constraint.is_inside(
                points=xyzs, buffer=state.radius
            )
            xyzs = xyzs[is_inside_mask]

        if state.bias:
            xyzs = state.bias(candidates=xyzs)

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
            f"Unable to find a starting point at {starting_xyz}"
            f"Unable to find a starting point at {starting_xyz}"
            "without overlapping particles. "
            "Check your `initial_point` argument."
        )
    elif state.volume_constraint:  # no initial point, but stay inside volume constraint
        xyzs = state.volume_constraint.sample_candidates(
            points=existing_points, n_candidates=100, buffer=state.radius + 0.1
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
            min_dist = -1*state.radius
        else:
            max_dist = np.max(existing_points, axis=0) 
            min_dist = np.min(existing_points, axis=0) # TODO: Update seed based on initial_point
        xyzs = state.rng.uniform(low=min_dist, high=max_dist, size=(100,3))
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
            assert 'low' in kwargs
            assert 'high' in kwargs
        elif distributionStr.lower() == "normal":
            self.sampler = self.rng.normal
            assert 'loc' in kwargs
            assert 'scale' in kwargs
        elif distributionStr == "choice":
            self.sampler = self.rng.choice
            assert 'a' in kwargs # p is not required
        else:
            raise NotImplementedError(f"Sample Distribution {distributionStr} not supported.")
        self.kwargs = kwargs

    def sample(self, size=None):
        return self.sampler(size=size, **self.kwargs)


def generate_trials(state):
    """Use normal or uniform sampling on angles, uniform sampling on radius."""
    thetas = state.angles.sample(
        size=state.trial_batch_size
    ).astype(np.float32)
    r = state.rng.uniform(-0.5, 0.5, size=(state.trial_batch_size, 3)).astype(
        np.float32
    )
    return thetas, r

def generate_crosslink_init_points():
    pass

