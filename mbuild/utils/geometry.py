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
