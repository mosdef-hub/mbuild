from numpy import *
from numpy.linalg import norm, svd, inv


__all__ = ['rotate_around_x', 'rotate_around_y', 'rotate_around_z',
           'equivalence_transform', 'translate', 'translate_to',
           'x_axis_transform', 'y_axis_transform', 'z_axis_transform']


class CoordinateTransform(object):
    """  """
    def __init__(self, T=None):
        if T is None:
            T = eye(4)

        self.T = T
        self.Tinv = inv(T)

    def apply_to(self, A):
        """Apply the coordinate transformation to points in A. """
        if A.ndim == 1:
            A = np.expand_dims(A, axis=0)
        rows, cols = A.shape
        A_new = hstack([A, ones((rows, 1))])

        A_new = transpose(self.T.dot(transpose(A_new)))
        return A_new[:, 0:cols]


class Translation(CoordinateTransform):
    """Cartesian translation. """
    def __init__(self, P):
        T = eye(4)
        T[0, 3] = P[0]
        T[1, 3] = P[1]
        T[2, 3] = P[2]
        super(Translation, self).__init__(T)


class RotationAroundZ(CoordinateTransform):
    """Rotation around the z-axis. """
    def __init__(self, theta):
        T = eye(4)
        T[0, 0] = cos(theta)
        T[0, 1] = -sin(theta)
        T[1, 0] = sin(theta)
        T[1, 1] = cos(theta)
        super(RotationAroundZ, self).__init__(T)


class RotationAroundY(CoordinateTransform):
    """Rotation around the y-axis. """
    def __init__(self, theta):
        T = eye(4)
        T[0, 0] = cos(theta)
        T[0, 2] = sin(theta)
        T[2, 0] = -sin(theta)
        T[2, 2] = cos(theta)
        super(RotationAroundY, self).__init__(T)


class RotationAroundX(CoordinateTransform):
    """Rotation around the x-axis. """
    def __init__(self, theta):
        T = eye(4)
        T[1, 1] = cos(theta)
        T[1, 2] = -sin(theta)
        T[2, 1] = sin(theta)
        T[2, 2] = cos(theta)
        super(RotationAroundX, self).__init__(T)


class Rotation(CoordinateTransform):
    """Rotation around vector by angle theta. """
    def __init__(self, theta, around):
        assert around.size == 3

        T = eye(4)

        s = sin(theta)
        c = cos(theta)
        t = 1 - c

        n = around / norm(around)

        x = n[0]
        y = n[1]
        z = n[2]
        m = array([
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c]])
        T[0:3, 0:3] = m
        super(Rotation, self).__init__(T)


class ChangeOfBasis(CoordinateTransform):
    """ """
    def __init__(self, basis, origin=None):
        assert shape(basis) == (3, 3)
        if origin is None:
            origin = array([0.0, 0.0, 0.0])

        T = eye(4)

        T[0:3, 0:3] = basis
        T = inv(T)

        T[0:3, 3:4] = -array([origin]).transpose()
        super(ChangeOfBasis, self).__init__(T)


class AxisTransform(CoordinateTransform):
    """ """
    def __init__(self, new_origin=None, point_on_x_axis=None,
                 point_on_xy_plane=None):
        if new_origin is None:
            new_origin = array([0.0, 0.0, 0.0])
        if point_on_x_axis is None:
            point_on_x_axis = array([1.0, 0.0, 0.0])
        if point_on_xy_plane is None:
            point_on_xy_plane = array([1.0, 1.0, 0.0])
        # Change the basis such that p1 is the origin, p2 is on the x axis and
        # p3 is in the xy plane.
        p1 = new_origin
        p2 = point_on_x_axis  # positive x axis
        p3 = point_on_xy_plane  # positive y part of the x axis

        # The direction vector of our new x axis.
        newx = unit_vector(p2 - p1)
        p3_u = unit_vector(p3 - p1)
        newz = unit_vector(cross(newx, p3_u))
        newy = cross(newz, newx)

        # Translation that moves new_origin to the origin.
        T_tr = eye(4)
        T_tr[0:3, 3:4] = -array([p1]).transpose()

        # Rotation that moves newx to the x axis, newy to the y axis, newz to
        # the z axis.
        B = eye(4)
        B[0:3, 0:3] = vstack((newx, newy, newz))

        # The concatentaion of translation and rotation.
        B_tr = dot(B, T_tr)

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
        rows, cols = shape(A)
        centroid_A = mean(A, axis=0)
        centroid_B = mean(B, axis=0)
        centroid_A.shape = (1, 3)
        centroid_B.shape = (1, 3)

        H = zeros((3, 3), dtype=float)

        for i in range(rows):
            H = H + transpose(A[i, :] - centroid_A).dot((B[i, :] - centroid_B))

        U, s, V = svd(H)
        V = transpose(V)
        R = V.dot(transpose(U))

        C_A = eye(3)
        C_A = vstack([hstack([C_A, transpose(centroid_A) * -1.0]),
                      array([0, 0, 0, 1])])

        R_new = eye(3)
        R_new = vstack([hstack([R, array([[0], [0], [0]])]),
                        array([0, 0, 0, 1])])

        C_B = eye(3)
        C_B = vstack([hstack([C_B, transpose(centroid_B)]),
                      array([0, 0, 0, 1])])

        T = C_B.dot(R_new).dot(C_A)

        super(RigidTransform, self).__init__(T)


def unit_vector(vector):
    """Returns the unit vector of the vector. """
    return vector / norm(vector)


def vec_angle(v1, v2):
    """Returns the angle in radians between vectors 'v1' and 'v2'

    >>> vec_angle((1, 0, 0), (0, 1, 0))
    1.5707963267948966
    >>> vec_angle((1, 0, 0), (1, 0, 0))
    0.0
    >>> vec_angle((1, 0, 0), (-1, 0, 0))
    3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    d = dot(v1_u, v2_u)
    if abs(d - 1.0) < 0.000000001:
        angle = 0.0
    else:
        angle = arccos(d)
    if isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return pi
    return angle


def _write_back_atom_positions(compound, arrnx3):
    from mbuild.atom import Atom
    if isinstance(compound, Atom):
        compound.pos = squeeze(arrnx3)
    else:
        for atom, coords in zip(compound.yield_atoms(), arrnx3):
            atom.pos = coords


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
    from mbuild.atom import Atom
    from mbuild.compound import Compound
    self_points = array([])
    self_points.shape = (0, 3)
    other_points = array([])
    other_points.shape = (0, 3)

    for pair in equiv:
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise Exception('Equivalence pair not a 2-tuple')
        if not ((isinstance(pair[0], Compound) and isinstance(pair[1], Compound)) or
                (isinstance(pair[0], Atom) and isinstance(pair[1], Atom))):
            # TODO: is the Atom and Atom comparison necessary?
            raise Exception('Equivalence pair type mismatch: pair[0] is a {0} '
                            'and pair[1] is a {1}'.format(type(pair[0]), type(pair[1])))

        if isinstance(pair[0], Atom):
            self_points = vstack([self_points, pair[0].pos])
            other_points = vstack([other_points, pair[1].pos])
        if isinstance(pair[0], Compound):
            for atom0 in pair[0].yield_atoms():
                self_points = vstack([self_points, atom0.pos])
            for atom1 in pair[1].yield_atoms():
                other_points = vstack([other_points, atom1.pos])

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
    from mbuild.port import Port
    from mbuild.bond import Bond
    if isinstance(from_positions, (list, tuple)) and isinstance(to_positions, (list, tuple)):
        equivalence_pairs = zip(from_positions, to_positions)
    elif isinstance(from_positions, Port) and isinstance(to_positions, Port):
        equivalence_pairs = _choose_correct_port(from_positions, to_positions)
    else:
        equivalence_pairs = [(from_positions, to_positions)]

    T = _create_equivalence_transform(equivalence_pairs)
    atom_positions = compound.xyz_with_ports
    atom_positions = T.apply_to(atom_positions)
    _write_back_atom_positions(compound, atom_positions)

    if add_bond:
        if isinstance(from_positions, Port) and isinstance(to_positions, Port):
            compound.add(Bond(from_positions, to_positions))


def _choose_correct_port(from_port, to_port):
    """Chooses the direction when using an equivalence transform on two Ports.

    Each Port object actually contains 2 sets of 4 atoms, either of which can be
    used to make a connection with an equivalence transform. This function
    chooses the set of 4 atoms that makes the anchor atoms not overlap which is
    the intended behavior for most use-cases.

    TODO: -Increase robustness for cases where the anchors are a different
           distance from their respective ports.
          -Provide options in `equivalence_transform` to override this behavior.

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
    T = _create_equivalence_transform([(from_port.up, to_port.up)])
    new_position = T.apply_to(array(from_port.anchor.pos, ndmin=2))

    dist_between_anchors_up_up = norm(new_position[0] - to_port.anchor.pos)

    # Then matching a 'down' with an 'up' port.
    T = _create_equivalence_transform([(from_port.down, to_port.up)])
    new_position = T.apply_to(array(from_port.anchor.pos, ndmin=2))

    # Determine which transform places the anchors further away from each other.
    dist_between_anchors_down_up = norm(new_position[0] - to_port.anchor.pos)
    difference_between_distances = dist_between_anchors_down_up - dist_between_anchors_up_up

    if difference_between_distances > 0:
        correct_port = from_port.down
    else:
        correct_port = from_port.up
    return [(correct_port, to_port.up)]


def translate(compound, v):
    atom_positions = compound.xyz_with_ports
    atom_positions = Translation(v).apply_to(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def translate_to(compound, v):
    atom_positions = compound.xyz_with_ports
    atom_positions -= compound.center
    atom_positions = Translation(v).apply_to(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def rotate_around_z(compound, theta):
    atom_positions = compound.xyz_with_ports
    atom_positions = RotationAroundZ(theta).apply_to(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def rotate_around_y(compound, theta):
    atom_positions = compound.xyz_with_ports
    atom_positions = RotationAroundY(theta).apply_to(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def rotate_around_x(compound, theta):
    atom_positions = compound.xyz_with_ports
    atom_positions = RotationAroundX(theta).apply_to(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def x_axis_transform(compound, new_origin=None,
                     point_on_x_axis=None,
                     point_on_xy_plane=None):
    """

    Parameters
    ----------
    compound:
    new_origin:
    point_on_x_axis:
    point_on_xy_plane:

    """
    from mbuild.compound import Atom

    if isinstance(new_origin, Atom):
        new_origin = new_origin.pos
    elif new_origin is None:
        new_origin = array([0, 0, 0])
    if isinstance(point_on_x_axis, Atom):
        point_on_x_axis = point_on_x_axis.pos
    elif point_on_x_axis is None:
        point_on_x_axis = array([1.0, 0.0, 0.0])
    if isinstance(point_on_xy_plane, Atom):
        point_on_xy_plane = point_on_xy_plane.pos
    elif point_on_xy_plane is None:
        point_on_xy_plane = array([1.0, 1.0, 0.0])

    atom_positions = compound.xyz_with_ports
    transform = AxisTransform(new_origin=new_origin,
                              point_on_x_axis=point_on_x_axis,
                              point_on_xy_plane=point_on_xy_plane)
    atom_positions = transform.apply_to(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def y_axis_transform(compound, new_origin=None,
                     point_on_y_axis=None,
                     point_on_xy_plane=None):
    x_axis_transform(compound, new_origin=new_origin,
                     point_on_x_axis=point_on_y_axis,
                     point_on_xy_plane=point_on_xy_plane)
    rotate_around_z(compound, pi / 2)


def z_axis_transform(compound, new_origin=None,
                     point_on_z_axis=None,
                     point_on_zx_plane=None):
    x_axis_transform(compound, new_origin=new_origin,
                     point_on_x_axis=point_on_z_axis,
                     point_on_xy_plane=point_on_zx_plane)
    rotate_around_y(compound, pi * 3 / 2)


def revolve_around_x(compound, theta):
    center_pos = compound.center
    translate(compound, -center_pos)
    rotate_around_x(compound, theta)
    translate(compound, center_pos)


def revolve_around_y(compound, theta):
    center_pos = compound.center
    translate(compound, -center_pos)
    rotate_around_y(compound, theta)
    translate(compound, center_pos)


def revolve_around_z(compound, theta):
    center_pos = compound.center
    translate(compound, -center_pos)
    rotate_around_z(compound, theta)
    translate(compound, center_pos)
