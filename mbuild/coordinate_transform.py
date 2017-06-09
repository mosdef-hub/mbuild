from warnings import warn, simplefilter
simplefilter('always', DeprecationWarning)

import numpy as np
from numpy.linalg import norm, svd, inv
from mbuild.utils.decorators import deprecated


__all__ = ['rotate', 'rotate_around_x', 'rotate_around_y', 'rotate_around_z',
           'spin', 'spin_x', 'spin_y', 'spin_z',
           'force_overlap', 'translate', 'translate_to',
           'x_axis_transform', 'y_axis_transform', 'z_axis_transform',

           # Deprecated
           'equivalence_transform']


def force_overlap(move_this, from_positions, to_positions, add_bond=True):
    """Computes an affine transformation that maps the from_positions to the
    respective to_positions, and applies this transformation to the compound.

    Parameters
    ----------
    move_this : mb.Compound
        The Compound to be moved.
    from_positions : np.ndarray, shape=(n, 3), dtype=float
        Original positions.
    to_positions : np.ndarray, shape=(n, 3), dtype=float
        New positions.
    add_bond : bool, optional, default=True
        If `from_positions` and `to_positions` are `Ports`, create a bond
        between the two anchor atoms.

    """
    from mbuild.port import Port
    T = None
    if isinstance(from_positions, (list, tuple)) and isinstance(to_positions, (list, tuple)):
        equivalence_pairs = zip(from_positions, to_positions)
    elif isinstance(from_positions, Port) and isinstance(to_positions, Port):
        equivalence_pairs, T = _choose_correct_port(from_positions, to_positions)
        from_positions.used = True
        to_positions.used = True
    else:
        equivalence_pairs = [(from_positions, to_positions)]

    if not T:
        T = _create_equivalence_transform(equivalence_pairs)
    atom_positions = move_this.xyz_with_ports
    atom_positions = T.apply_to(atom_positions)
    move_this.xyz_with_ports = atom_positions

    if add_bond:
        if isinstance(from_positions, Port) and isinstance(to_positions, Port):
            if not from_positions.anchor or not to_positions.anchor:
                warn("Attempting to form bond from port that has no anchor")
            else:
                from_positions.anchor.parent.add_bond((from_positions.anchor, to_positions.anchor))
                to_positions.anchor.parent.add_bond((from_positions.anchor, to_positions.anchor))


class CoordinateTransform(object):
    """  """
    def __init__(self, T=None):
        if T is None:
            T = np.eye(4)

        self.T = T
        self.Tinv = inv(T)

    def apply_to(self, A):
        """Apply the coordinate transformation to points in A. """
        if A.ndim == 1:
            A = np.expand_dims(A, axis=0)
        rows, cols = A.shape
        A_new = np.hstack([A, np.ones((rows, 1))])

        A_new = np.transpose(self.T.dot(np.transpose(A_new)))
        return A_new[:, 0:cols]


class Translation(CoordinateTransform):
    """Cartesian translation. """
    def __init__(self, P):
        T = np.eye(4)
        T[0, 3] = P[0]
        T[1, 3] = P[1]
        T[2, 3] = P[2]
        super(Translation, self).__init__(T)


class RotationAroundZ(CoordinateTransform):
    """Rotation around the z-axis. """
    def __init__(self, theta):
        T = np.eye(4)
        T[0, 0] = np.cos(theta)
        T[0, 1] = -np.sin(theta)
        T[1, 0] = np.sin(theta)
        T[1, 1] = np.cos(theta)
        super(RotationAroundZ, self).__init__(T)


class RotationAroundY(CoordinateTransform):
    """Rotation around the y-axis. """
    def __init__(self, theta):
        T = np.eye(4)
        T[0, 0] = np.cos(theta)
        T[0, 2] = np.sin(theta)
        T[2, 0] = -np.sin(theta)
        T[2, 2] = np.cos(theta)
        super(RotationAroundY, self).__init__(T)


class RotationAroundX(CoordinateTransform):
    """Rotation around the x-axis. """
    def __init__(self, theta):
        T = np.eye(4)
        T[1, 1] = np.cos(theta)
        T[1, 2] = -np.sin(theta)
        T[2, 1] = np.sin(theta)
        T[2, 2] = np.cos(theta)
        super(RotationAroundX, self).__init__(T)


class Rotation(CoordinateTransform):
    """Rotation around vector by angle theta. """
    def __init__(self, theta, around):
        assert around.size == 3

        T = np.eye(4)

        s = np.sin(theta)
        c = np.cos(theta)
        t = 1 - c

        n = around / norm(around)

        x = n[0]
        y = n[1]
        z = n[2]
        m = np.array([
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c]])
        T[0:3, 0:3] = m
        super(Rotation, self).__init__(T)


class ChangeOfBasis(CoordinateTransform):
    """ """
    def __init__(self, basis, origin=None):
        assert np.shape(basis) == (3, 3)
        if origin is None:
            origin = np.array([0.0, 0.0, 0.0])

        T = np.eye(4)

        T[0:3, 0:3] = basis
        T = inv(T)

        T[0:3, 3:4] = -np.array([origin]).transpose()
        super(ChangeOfBasis, self).__init__(T)


class AxisTransform(CoordinateTransform):
    """ """
    def __init__(self, new_origin=None, point_on_x_axis=None,
                 point_on_xy_plane=None):
        if new_origin is None:
            new_origin = np.array([0.0, 0.0, 0.0])
        if point_on_x_axis is None:
            point_on_x_axis = np.array([1.0, 0.0, 0.0])
        if point_on_xy_plane is None:
            point_on_xy_plane = np.array([1.0, 1.0, 0.0])
        # Change the basis such that p1 is the origin, p2 is on the x axis and
        # p3 is in the xy plane.
        p1 = new_origin
        p2 = point_on_x_axis  # positive x axis
        p3 = point_on_xy_plane  # positive y part of the x axis

        # The direction vector of our new x axis.
        newx = unit_vector(p2 - p1)
        p3_u = unit_vector(p3 - p1)
        newz = unit_vector(np.cross(newx, p3_u))
        newy = np.cross(newz, newx)

        # Translation that moves new_origin to the origin.
        T_tr = np.eye(4)
        T_tr[0:3, 3:4] = -np.array([p1]).transpose()

        # Rotation that moves newx to the x axis, newy to the y axis, newz to
        # the z axis.
        B = np.eye(4)
        B[0:3, 0:3] = np.vstack((newx, newy, newz))

        # The concatentaion of translation and rotation.
        B_tr = np.dot(B, T_tr)

        super(AxisTransform, self).__init__(B_tr)


class RigidTransform(CoordinateTransform):
    """Computes the rigid transformation that maps points A to points B.

    See http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.173.2196&rep=rep1&type=pdf
    """

    def __init__(self, A, B):
        """

        Parameters
        ----------
        A : np.ndarray, shape=(n, 3), dtype=float
            Points in source coordinate system.
        B : np.ndarray, shape=(n, 3), dtype=float
            Points in destination coordinate system.

        """
        rows, _ = np.shape(A)
        centroid_A = np.mean(A, axis=0)
        centroid_B = np.mean(B, axis=0)
        centroid_A.shape = (1, 3)
        centroid_B.shape = (1, 3)

        H = np.zeros((3, 3), dtype=float)

        for i in range(rows):
            H = H + np.transpose(A[i, :] - centroid_A).dot((B[i, :] - centroid_B))

        U, _, V = svd(H)
        V = np.transpose(V)
        R = V.dot(np.transpose(U))

        C_A = np.eye(3)
        C_A = np.vstack([np.hstack([C_A, np.transpose(centroid_A) * -1.0]),
                         np.array([0, 0, 0, 1])])

        R_new = np.vstack([np.hstack([R, np.array([[0], [0], [0]])]),
                           np.array([0, 0, 0, 1])])

        C_B = np.eye(3)
        C_B = np.vstack([np.hstack([C_B, np.transpose(centroid_B)]),
                        np.array([0, 0, 0, 1])])

        T = C_B.dot(R_new).dot(C_A)

        super(RigidTransform, self).__init__(T)


def unit_vector(v):
    """Returns the unit vector of the vector. """
    return v / norm(v)


def angle(u, v, w=None):
    """Returns the angle in radians between two vectors. """
    if w != None:
        u = u - v
        v = w - v
    c = np.dot(u, v) / norm(u) / norm(v)
    return np.arccos(np.clip(c, -1, 1))


def _create_equivalence_transform(equiv):
    """Compute an equivalence transformation that transforms this compound
    to another compound's coordinate system.

    Parameters
    ----------
    equiv : np.ndarray, shape=(n, 3), dtype=float
        Array of equivalent points.

    Returns
    -------
    T : CoordinateTransform
        Transform that maps this point cloud to the other point cloud's
        coordinates system.

    """
    from mbuild.compound import Compound
    self_points = np.array([])
    self_points.shape = (0, 3)
    other_points = np.array([])
    other_points.shape = (0, 3)

    for pair in equiv:
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise ValueError('Equivalence pair not a 2-tuple')
        if not (isinstance(pair[0], Compound) and isinstance(pair[1], Compound)):
            raise ValueError('Equivalence pair type mismatch: pair[0] is a {0} '
                             'and pair[1] is a {1}'.format(type(pair[0]),
                                                           type(pair[1])))

        # TODO: vstack is slow, replace with list concatenation
        if not pair[0].children:
            self_points = np.vstack([self_points, pair[0].pos])
            other_points = np.vstack([other_points, pair[1].pos])
        else:
            for atom0 in pair[0]._particles(include_ports=True):
                self_points = np.vstack([self_points, atom0.pos])
            for atom1 in pair[1]._particles(include_ports=True):
                other_points = np.vstack([other_points, atom1.pos])
    T = RigidTransform(self_points, other_points)
    return T


def equivalence_transform(compound, from_positions, to_positions, add_bond=True):
    """Computes an affine transformation that maps the from_positions to the
    respective to_positions, and applies this transformation to the compound.

    Parameters
    ----------
    compound : mb.Compound
        The Compound to be transformed.
    from_positions : np.ndarray, shape=(n, 3), dtype=float
        Original positions.
    to_positions : np.ndarray, shape=(n, 3), dtype=float
        New positions.

    """
    warn('The `equivalence_transform` function is being phased out in favor of'
         ' `force_overlap`.', DeprecationWarning)
    from mbuild.port import Port
    T = None
    if isinstance(from_positions, (list, tuple)) and isinstance(to_positions, (list, tuple)):
        equivalence_pairs = zip(from_positions, to_positions)
    elif isinstance(from_positions, Port) and isinstance(to_positions, Port):
        equivalence_pairs, T = _choose_correct_port(from_positions, to_positions)
        from_positions.used = True
        to_positions.used = True
    else:
        equivalence_pairs = [(from_positions, to_positions)]

    if not T:
        T = _create_equivalence_transform(equivalence_pairs)
    atom_positions = compound.xyz_with_ports
    atom_positions = T.apply_to(atom_positions)
    compound.xyz_with_ports = atom_positions

    if add_bond:
        if isinstance(from_positions, Port) and isinstance(to_positions, Port):
            if not from_positions.anchor or not to_positions.anchor:
                # TODO: I think warnings is undefined here
                warn("Attempting to form bond from port that has no anchor")
            else:
                from_positions.anchor.parent.add_bond((from_positions.anchor, to_positions.anchor))
                to_positions.anchor.parent.add_bond((from_positions.anchor, to_positions.anchor))


def _choose_correct_port(from_port, to_port):
    """Chooses the direction when using an equivalence transform on two Ports.

    Each Port object actually contains 2 sets of 4 atoms, either of which can be
    used to make a connection with an equivalence transform. This function
    chooses the set of 4 atoms that makes the anchor atoms not overlap which is
    the intended behavior for most use-cases.

    TODO: -Increase robustness for cases where the anchors are a different
           distance from their respective ports.
          -Provide options in `force_overlap` to override this behavior.

    Parameters
    ----------
    from_port : mb.Port
    to_port : mb.Port

    Returns
    -------
    equivalence_pairs : tuple of Ports, shape=(2,)
        Technically, a tuple of the Ports' sub-Compounds ('up' or 'down')
        that are used to make the correct connection between components.

    """
    # First we try matching the two 'up' ports.
    T1 = _create_equivalence_transform([(from_port['up'], to_port['up'])])
    new_position = T1.apply_to(np.array(from_port.anchor.pos, ndmin=2))

    dist_between_anchors_up_up = norm(new_position[0] - to_port.anchor.pos)

    # Then matching a 'down' with an 'up' port.
    T2 = _create_equivalence_transform([(from_port['down'], to_port['up'])])
    new_position = T2.apply_to(np.array(from_port.anchor.pos, ndmin=2))

    # Determine which transform places the anchors further away from each other.
    dist_between_anchors_down_up = norm(new_position[0] - to_port.anchor.pos)
    difference_between_distances = dist_between_anchors_down_up - dist_between_anchors_up_up

    if difference_between_distances > 0:
        correct_port = from_port['down']
        T = T2
    else:
        correct_port = from_port['up']
        T = T1
    return [(correct_port, to_port['up'])], T


warning_message = 'Please use Compound.translate()'
@deprecated(warning_message)
def translate(compound, pos):
    """Translate a compound by a vector.

    Parameters
    ----------
    compound : mb.Compound
        The compound being translated.
    pos : np.ndarray, shape=(3,), dtype=float
        The vector to translate the compound by.

    """
    atom_positions = compound.xyz_with_ports
    atom_positions = Translation(pos).apply_to(atom_positions)
    compound.xyz_with_ports = atom_positions

warning_message = 'Please use Compound.translate_to()'
@deprecated(warning_message)
def translate_to(compound, pos):
    """Translate a compound to a coordinate.

    Parameters
    ----------
    compound : mb.Compound
        The compound being translated.
    pos : np.ndarray, shape=(3,), dtype=float
        The coordinate to translate the compound to.

    """
    atom_positions = compound.xyz_with_ports
    atom_positions -= compound.center
    atom_positions = Translation(pos).apply_to(atom_positions)
    compound.xyz_with_ports = atom_positions


def _translate(coordinates, by):
    """Translate a set of coordinates by a vector.

    Parameters
    ----------
    coordinates : np.ndarray, shape=(n,3), dtype=float
        The coordinates being translated.
    by : np.ndarray, shape=(3,), dtype=float
        The vector to translate the coordinates by.

    """
    return Translation(by).apply_to(coordinates)


def _translate_to(coordinates, to):
    """Translate a set of coordinates to a location.

    Parameters
    ----------
    coordinates : np.ndarray, shape=(n,3), dtype=float
        The coordinates being translated.
    to : np.ndarray, shape=(3,), dtype=float
        The new average position of the coordinates.

    """
    coordinates -= np.mean(coordinates, axis=0)
    return Translation(to).apply_to(coordinates)


def _rotate(coordinates, theta, around):
    """Rotate a set of coordinates around an arbitrary vector.

    Parameters
    ----------
    coordinates : np.ndarray, shape=(n,3), dtype=float
        The coordinates being rotated.
    theta : float
        The angle by which to rotate the coordinates, in radians.
    around : np.ndarray, shape=(3,), dtype=float
        The vector about which to rotate the coordinates.

    """
    around = np.asarray(around).reshape(3)
    if np.array_equal(around, np.zeros(3)):
        raise ValueError('Cannot rotate around a zero vector')
    return Rotation(theta, around).apply_to(coordinates)


warning_message = 'Please use Compound.rotate()'
@deprecated(warning_message)
def rotate(compound, theta, around):
    """Rotate a compound around an arbitrary vector.

    Parameters
    ----------
    compound : mb.Compound
        The compound being rotated.
    theta : float
        The angle by which to rotate the compound, in radians.
    around : np.ndarray, shape=(3,), dtype=float
        The vector about which to rotate the compound.

    """
    around = np.asarray(around).reshape(3)
    if np.array_equal(around, np.zeros(3)):
        raise ValueError('Cannot rotate around a zero vector')
    atom_positions = compound.xyz_with_ports
    atom_positions = Rotation(theta, around).apply_to(atom_positions)
    compound.xyz_with_ports = atom_positions


warning_message = 'Please use rotate(compound, theta, around=np.asarray([1, 0, 0]))'
@deprecated(warning_message)
def rotate_around_x(compound, theta):
    """Rotate a compound around the x axis.

    Parameters
    ----------
    compound : mb.Compound
        The compound being rotated.
    theta : float
        The angle by which to rotate the compound.

    """
    rotate(compound, theta, np.asarray([1, 0, 0]))


warning_message = 'Please use rotate(compound, theta, around=np.asarray([0, 1, 0]))'
@deprecated(warning_message)
def rotate_around_y(compound, theta):
    """Rotate a compound around the y axis.

    Parameters
    ----------
    compound : mb.Compound
        The compound being rotated.
    theta : float
        The angle by which to rotate the compound.

    """
    rotate(compound, theta, np.asarray([0, 1, 0]))


warning_message = 'Please use rotate(compound, theta, around=np.asarray([0, 0, 1]))'
@deprecated(warning_message)
def rotate_around_z(compound, theta):
    """Rotate a compound around the z axis.

    Parameters
    ----------
    compound : mb.Compound
        The compound being rotated.
    theta : float
        The angle by which to rotate the compound.

    """
    rotate(compound, theta, np.asarray([0, 0, 1]))


warning_message = 'Please use Compound.spin()'
@deprecated(warning_message)
def spin(compound, theta, around):
    """Rotate a compound in place around an arbitrary vector.

    Parameters
    ----------
    compound : mb.Compound
        The compound being rotated.
    theta : float
        The angle by which to rotate the compound, in radians.
    around : np.ndarray, shape=(3,), dtype=float
        The axis about which to spin the compound.

    """
    around = np.asarray(around).reshape(3)
    if np.array_equal(around, np.zeros(3)):
        raise ValueError('Cannot spin around a zero vector')
    center_pos = compound.center
    translate(compound, -center_pos)
    rotate(compound, theta, around)
    translate(compound, center_pos)


def _spin(coordinates, theta, around):
    """Rotate a set of coordinates in place around an arbitrary vector.

    Parameters
    ----------
    coordinates : np.ndarray, shape=(n,3), dtype=float
        The coordinates being spun.
    theta : float
        The angle by which to spin the coordinates, in radians.
    around : np.ndarray, shape=(3,), dtype=float
        The axis about which to spin the coordinates.

    """
    around = np.asarray(around).reshape(3)
    if np.array_equal(around, np.zeros(3)):
        raise ValueError('Cannot spin around a zero vector')
    center_pos = np.mean(coordinates, axis=0)
    coordinates -= center_pos
    coordinates = _rotate(coordinates, theta, around)
    coordinates += center_pos
    return coordinates


warning_message = 'Please use spin(compound, theta, around=np.asarray([1, 0, 0]))'
@deprecated(warning_message)
def spin_x(compound, theta):
    """Rotate a compound in place around the x axis.

    Parameters
    ----------
    compound : mb.Compound
        The compound being rotated.
    theta : float
        The angle by which to rotate the compound.

    """
    spin(compound, theta, np.asarray([1.0, 0.0, 0.0]))


warning_message = 'Please use spin(compound, theta, around=np.asarray([0, 1, 0]))'
@deprecated(warning_message)
def spin_y(compound, theta):
    """Rotate a compound in place around the y axis.

    Parameters
    ----------
    compound : mb.Compound
        The compound being rotated.
    theta : float
        The angle by which to rotate the compound.

    """
    spin(compound, theta, np.asarray([0.0, 1.0, 0.0]))


warning_message = 'Please use spin(compound, theta, around=np.asarray([0, 0, 1]))'
@deprecated(warning_message)
def spin_z(compound, theta):
    """Rotate a compound in place around the z axis.

    Parameters
    ----------
    compound : mb.Compound
        The compound being rotated.
    theta : float
        The angle by which to rotate the compound.

    """
    spin(compound, theta, np.asarray([0.0, 0.0, 1.0]))


def x_axis_transform(compound, new_origin=None,
                     point_on_x_axis=None,
                     point_on_xy_plane=None):
    """Move a compound such that the x-axis lies on specified points.

    Parameters
    ----------
    compound : mb.Compound
        The compound to move.
    new_origin : mb.Compound or list-like, optional, default=[0.0, 0.0, 0.0]
        Where to place the new origin of the coordinate system.
    point_on_x_axis : mb.Compound or list-like, optional, default=[1.0, 0.0, 0.0]
        A point on the new x-axis.
    point_on_xy_plane : mb.Compound, or list-like, optional, default=[1.0, 0.0, 0.0]
        A point on the new xy-plane.

    """
    import mbuild as mb

    if new_origin is None:
        new_origin = np.array([0, 0, 0])
    elif isinstance(new_origin, mb.Compound):
        new_origin = new_origin.pos
    elif isinstance(new_origin, (tuple, list,np.ndarray)):
        new_origin = np.asarray(new_origin)
    else:
        raise TypeError('x_axis_transform, y_axis_transform, and z_axis_transform only accept'
                        ' np.ndarrays, mb.Compounds, lists, or None of size 3 for the new_origin'
                        ' parameter. User passed type: {}.'.format(type(new_origin)))
    if point_on_x_axis is None:
        point_on_x_axis = np.array([1.0, 0.0, 0.0])
    elif isinstance(point_on_x_axis, mb.Compound):
        point_on_x_axis = point_on_x_axis.pos
    elif isinstance(point_on_x_axis, (list, tuple, np.ndarray)):
        point_on_x_axis = np.asarray(point_on_x_axis)
    else:
        raise TypeError('x_axis_transform, y_axis_transform, and z_axis_transform only accept'
                         ' np.ndarrays, mb.Compounds, lists, or None of size 3 for the'
                         ' point_on_x_axis parameter. User passed type: {}.'.format(type(point_on_x_axis)))    
    if point_on_xy_plane is None:
        point_on_xy_plane = np.array([1.0, 1.0, 0.0])
    elif isinstance(point_on_xy_plane, mb.Compound):
        point_on_xy_plane = point_on_xy_plane.pos
    elif isinstance(point_on_xy_plane, (list, tuple, np.ndarray)):
        point_on_xy_plane = np.asarray(point_on_xy_plane)
    else:
        raise TypeError('x_axis_transform, y_axis_transform, and z_axis_transform only accept'
                          ' np.ndarrays, mb.Compounds, lists, or None of size 3 for the'
                          ' point_on_xy_plane parameter. User passed type: {}.'.format(type(point_on_xy_plane)))

    atom_positions = compound.xyz_with_ports
    transform = AxisTransform(new_origin=new_origin,
                              point_on_x_axis=point_on_x_axis,
                              point_on_xy_plane=point_on_xy_plane)
    atom_positions = transform.apply_to(atom_positions)
    compound.xyz_with_ports = atom_positions


def y_axis_transform(compound, new_origin=None,
                     point_on_y_axis=None,
                     point_on_xy_plane=None):
    """Move a compound such that the y-axis lies on specified points.

    Parameters
    ----------
    compound : mb.Compound
        The compound to move.
    new_origin : mb.Compound or np.ndarray, optional, default=[0.0, 0.0, 0.0]
        Where to place the new origin of the coordinate system.
    point_on_y_axis : mb.Compound or np.ndarray, optional, default=[0.0, 1.0, 0.0]
        A point on the new x-axis.
    point_on_xy_plane : mb.Compound or np.ndarray, optional, default=[0.0, 1.0, 0.0]
        A point on the new xy-plane.

    """
    x_axis_transform(compound, new_origin=new_origin,
                     point_on_x_axis=point_on_y_axis,
                     point_on_xy_plane=point_on_xy_plane)
    rotate_around_z(compound, np.pi / 2)


def z_axis_transform(compound, new_origin=None,
                     point_on_z_axis=None,
                     point_on_zx_plane=None):
    """Move a compound such that the z-axis lies on specified points.

    Parameters
    ----------
    compound : mb.Compound
        The compound to move.
    new_origin : mb.Compound or np.ndarray, optional, default=[0.0, 0.0, 0.0]
        Where to place the new origin of the coordinate system.
    point_on_y_axis : mb.Compound or np.ndarray, optional, default=[0.0, 0.0, 1.0]
        A point on the new z-axis.
    point_on_xy_plane : mb.Compound or np.ndarray, optional, default=[0.0, 0.0, 1.0]
        A point on the new xz-plane.

    """
    x_axis_transform(compound, new_origin=new_origin,
                     point_on_x_axis=point_on_z_axis,
                     point_on_xy_plane=point_on_zx_plane)
    rotate_around_y(compound, np.pi * 3 / 2)
