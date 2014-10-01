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


def compute_inertia_tensor(traj, masses=None):
    """Calculates the moment of inertia tensor for each frame.

    Function adapted from MDAnalysis (https://code.google.com/p/mdanalysis/)

    Args:
        traj (Trajectory): Trajectory to compute center of mass for.
        masses (np.ndarray, optional):

    Returns:
        I (np.ndarray, shape=(n_frames, 3, 3): Moment of inertia tensor for each
        frame.
    """

    if not masses:
        try:
            masses = np.array([a.element.mass for a in traj.top.atoms])
        except:  # TODO: figure out types of errors to catch
            masses = np.ones(traj.n_atoms)
    else:
        assert masses.shape[0] == traj.xyz[0].shape[0]

    com = compute_center_of_mass(traj, masses)
    I = np.zeros(shape=(traj.n_frames, 3, 3))
    for n, frame in enumerate(traj.xyz):
        recenteredpos = frame - com[n]
        values = zip(masses, recenteredpos)
        Ixx = reduce(lambda t,a: t+a[0]*(a[1][1]*a[1][1]+a[1][2]*a[1][2]), values, 0.)
        Iyy = reduce(lambda t,a: t+a[0]*(a[1][0]*a[1][0]+a[1][2]*a[1][2]), values, 0.)
        Izz = reduce(lambda t,a: t+a[0]*(a[1][0]*a[1][0]+a[1][1]*a[1][1]), values, 0.)
        Ixy = Iyx = -1*reduce(lambda t,a: t+a[0]*a[1][0]*a[1][1], values, 0.)
        Ixz = Izx = -1*reduce(lambda t,a: t+a[0]*a[1][0]*a[1][2], values, 0.)
        Iyz = Izy = -1*reduce(lambda t,a: t+a[0]*a[1][1]*a[1][2], values, 0.)
        I[n] = np.array([[Ixx, Ixy, Ixz], [Iyx, Iyy, Iyz], [Izx, Izy, Izz]])
    return I