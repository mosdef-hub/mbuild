"""mBuild utilites for geometrical operations."""

import numpy as np

import mbuild as mb
from mbuild.coordinate_transform import angle


def calc_dihedral(point1, point2, point3, point4):
    """Calculate a dihedral angle.

    Here, two planes are defined by (point1, point2, point3) and (point2,
    point3, point4). The angle between them is returned.

    Parameters
    ----------
    point1, point2, point3, point4 : array-like, shape=(3,), dtype=float
        Four points that define two planes

    Returns
    -------
    float
        The dihedral angle between the two planes defined by the four points.
    """
    points = np.array([point1, point2, point3, point4])
    x = np.cross(points[1] - points[0], points[2] - points[1])
    y = np.cross(points[2] - points[1], points[3] - points[2])
    return angle(x, y)


def coord_shift(xyz, box):
    """Ensure that coordinates are -L/2, L/2.

    Checks if coordinates are -L/2, L/2 and then shifts coordinates if
    necessary. For example, if coordinates are 0, L, then a shift is applied to
    move coordinates to -L/2, L/2. If a shift is not necessary, the points are
    returned unmodified.

    Parameters
    ----------
    xyz : numpy.array, shape=(N,3)
        Coordinates
    box : numpy.array, shape=(N,3)
        Array specifing the box lengths, e.g., [Lx, Ly, Lz]

    Returns
    -------
    xyz : numpy.array, shape=(N,3)
        Shifted coordinates
    """
    box = np.asarray(box)
    assert box.shape == (3,)

    box_max = box / 2.0
    box_min = -box_max
    # Shift all atoms
    if np.greater(xyz, box_max).any():
        xyz -= box_max
    elif np.less(xyz, box_min).any():
        xyz += box_max

    return xyz


def wrap_coords(xyz, box, mins=None):
    """Wrap coordinates inside box.

    Parameters
    ----------
    xyz : numpy.array, shape=(N,3),
        Coordinates
    box : numpy.array or list or mb.Box
        Array or list should have shape (3,) corresponding to box lengths.
        If array or list is passed, box is assumed to be positive octant
        If mb.box is passed, box can be arbitrarily centered

    Returns
    -------
    wrap_xyz : numpy.array, shape=(N,3)
        wrapped coordinates

    Notes
    -----
    Currently only supports orthorhombic boxes
    """
    if not isinstance(box, mb.Box):
        box_arr = np.asarray(box)
        assert box_arr.shape == (3,)

        wrap_xyz = xyz - 1 * np.floor_divide(xyz, box_arr) * box_arr
    else:
        xyz = xyz - mins
        wrap_xyz = (
            xyz
            - (np.floor_divide(xyz, np.asarray(box.lengths)) * np.asarray(box.lengths))
            + mins
        )

    return wrap_xyz
