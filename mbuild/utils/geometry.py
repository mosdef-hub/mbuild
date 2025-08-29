"""mBuild utilites for geometrical operations."""

import numpy as np

from mbuild.box import Box
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


def bounding_box(xyz):
    """Find the bounding box from a set of coordinates."""
    # case where only 1 particle exists
    is_one_particle = False
    if xyz.shape[0] == 1:
        is_one_particle = True

    # are any columns all equalivalent values?
    # an example of this would be a planar molecule
    # example: all z values are 0.0
    # from: https://stackoverflow.com/a/14860884
    # steps: create mask array comparing first value in each column
    # use np.all with axis=0 to do row columnar comparision
    has_dimension = [True, True, True]
    if not is_one_particle:
        missing_dimensions = np.all(
            np.isclose(xyz, xyz[0, :], atol=1e-2),
            axis=0,
        )
        for i, truthy in enumerate(missing_dimensions):
            has_dimension[i] = not truthy

    if is_one_particle:
        v1 = np.asarray([[1.0, 0.0, 0.0]])
        v2 = np.asarray([[0.0, 1.0, 0.0]])
        v3 = np.asarray([[0.0, 0.0, 1.0]])
    else:
        maxs = xyz.max(axis=0)
        mins = xyz.min(axis=0)
        v1 = np.asarray((maxs[0] - mins[0], 0.0, 0.0))
        v2 = np.asarray((0.0, maxs[1] - mins[1], 0.0))
        v3 = np.asarray((0.0, 0.0, maxs[2] - mins[2]))
    vecs = [v1, v2, v3]

    # handle any missing dimensions (planar molecules)
    for i, dim in enumerate(has_dimension):
        if not dim:
            vecs[i][i] = 0.1
    return vecs


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
    if not isinstance(box, Box):
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
