from __future__ import division

from functools import reduce

import numpy as np


def compute_center_of_mass(traj, masses=None):
    """Compute the center of mass for each frame.

    Note: This function alters the equivalent mdtraj by adding the option to
    provide an array of masses.

    Args:
        traj (Trajectory): Trajectory to compute center of mass for.
        masses (np.ndarray, optional):

    Returns:
        com (np.ndarray, shape=(n_frames, 3)): Coordinates of the center of mass
        for each frame.

    """

    com = np.zeros((traj.n_frames, 3))
    if masses is None:
        try:
            masses = np.array([a.element.mass for a in traj.top.atoms])
        except:
            masses = np.ones(traj.n_atoms)
    else:
        assert masses.shape[0] == traj.xyz[0].shape[0]

    masses /= masses.sum()

    for i, x in enumerate(traj.xyz):
        com[i, :] = x.astype('float64').T.dot(masses)
    return com

