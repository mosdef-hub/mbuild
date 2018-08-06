import numpy as np

from mbuild.coordinate_transform import angle

def calc_dihedral(point1, point2, point3, point4):
    """Calculates a dihedral angle

    Here, two planes are defined by (point1, point2, point3) and
    (point2, point3, point4). The angle between them is returned.

    Parameters
    ----------
    point1, point2, point3, point4 : array-like, shape=(3,), dtype=float
        Four points that define two planes

    Returns
    -------
    float
        The dihedral angle between the two planes defined by the four
        points.
    """
    points = np.array([point1, point2, point3, point4])
    x = np.cross(points[1] - points[0], points[2] - points[1])
    y = np.cross(points[2] - points[1], points[3] - points[2])
    return angle(x, y)


def shift_coords(xyz, box):
    """Ensures that coordinates are -L/2, L/2

    Checks if coordinates are -L/2, L/2 and then shifts coordinates
    if necessary. For example, if coordinates are 0, L, then a shift
    is applied to move coordinates to -L/2, L/2.

    Parameters
    ----------
    xyz : numpy.array of points with shape N x 3
    box : numpy.array specifing the size of box ie [Lx, Ly, Lz]
    """

    # Check if we need to shift things to -L/2, L/2
    box_max = structure.box[:3]/2.
    wrap = np.greater(xyz, box_max).any()
    # Shift all atoms by box_min
    if wrap:
        xyz -= box_max

