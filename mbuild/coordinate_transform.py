"""Coordinate transformation functions."""

from warnings import simplefilter, warn

simplefilter("always", DeprecationWarning)

import numpy as np
from numpy.linalg import inv, norm, svd

__all__ = [
    "force_overlap",
    "x_axis_transform",
    "y_axis_transform",
    "z_axis_transform",
    "equivalence_transform",
]


def force_overlap(
    move_this, from_positions, to_positions, add_bond=True, reset_labels=False
):
    """Move a Compound such that a position overlaps with another.

    Computes an affine transformation that maps the from_positions to the
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
    reset_labels : bool, optional, default=False
        If True, the Compound labels will be reset, renumbered using the Compound.reset_labels methods
    """
    from mbuild.port import Port

    T = None
    if isinstance(from_positions, (list, tuple)) and isinstance(
        to_positions, (list, tuple)
    ):
        equivalence_pairs = zip(from_positions, to_positions)
    elif isinstance(from_positions, Port) and isinstance(to_positions, Port):
        equivalence_pairs, T = _choose_correct_port(
            from_positions, to_positions
        )
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
                from_positions.anchor.parent.add_bond(
                    (from_positions.anchor, to_positions.anchor)
                )
                to_positions.anchor.parent.add_bond(
                    (from_positions.anchor, to_positions.anchor)
                )
                from_positions.anchor.parent.remove(
                    from_positions, reset_labels=reset_labels
                )
                to_positions.anchor.parent.remove(
                    to_positions, reset_labels=reset_labels
                )


class CoordinateTransform(object):
    """Coordinate transforms."""

    def __init__(self, T=None):
        if T is None:
            T = np.eye(4)

        self.T = T
        self.Tinv = inv(T)

    def apply_to(self, A):
        """Apply the coordinate transformation to points in A."""
        if A.ndim == 1:
            A = np.expand_dims(A, axis=0)
        rows, cols = A.shape
        A_new = np.hstack([A, np.ones((rows, 1))])

        A_new = np.transpose(self.T.dot(np.transpose(A_new)))
        return A_new[:, 0:cols]


class Translation(CoordinateTransform):
    """Cartesian translation."""

    def __init__(self, P):
        T = np.eye(4)
        T[0, 3] = P[0]
        T[1, 3] = P[1]
        T[2, 3] = P[2]
        super(Translation, self).__init__(T)


class RotationAroundZ(CoordinateTransform):
    """Rotation around the z-axis."""

    def __init__(self, theta):
        T = np.eye(4)
        T[0, 0] = np.cos(theta)
        T[0, 1] = -np.sin(theta)
        T[1, 0] = np.sin(theta)
        T[1, 1] = np.cos(theta)
        super(RotationAroundZ, self).__init__(T)


class RotationAroundY(CoordinateTransform):
    """Rotation around the y-axis."""

    def __init__(self, theta):
        T = np.eye(4)
        T[0, 0] = np.cos(theta)
        T[0, 2] = np.sin(theta)
        T[2, 0] = -np.sin(theta)
        T[2, 2] = np.cos(theta)
        super(RotationAroundY, self).__init__(T)


class RotationAroundX(CoordinateTransform):
    """Rotation around the x-axis."""

    def __init__(self, theta):
        T = np.eye(4)
        T[1, 1] = np.cos(theta)
        T[1, 2] = -np.sin(theta)
        T[2, 1] = np.sin(theta)
        T[2, 2] = np.cos(theta)
        super(RotationAroundX, self).__init__(T)


class Rotation(CoordinateTransform):
    """Rotation around vector by angle theta."""

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
        m = np.array(
            [
                [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
                [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
                [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
            ]
        )
        T[0:3, 0:3] = m
        super(Rotation, self).__init__(T)


class ChangeOfBasis(CoordinateTransform):
    """Convert the basis of coordinates to another basis."""

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
    """Axis transform."""

    def __init__(
        self, new_origin=None, point_on_x_axis=None, point_on_xy_plane=None
    ):
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

    See http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.173.2196&rep=
    rep1&type=pdf

    Parameters
    ----------
    A : np.ndarray, shape=(n, 3), dtype=float
        Points in source coordinate system.
    B : np.ndarray, shape=(n, 3), dtype=float
        Points in destination coordinate system.
    """

    def __init__(self, A, B):
        rows, _ = np.shape(A)
        centroid_A = np.mean(A, axis=0)
        centroid_B = np.mean(B, axis=0)
        centroid_A.shape = (1, 3)
        centroid_B.shape = (1, 3)

        H = np.zeros((3, 3), dtype=float)

        for i in range(rows):
            H = H + np.transpose(A[i, :] - centroid_A).dot(
                (B[i, :] - centroid_B)
            )

        U, _, V = svd(H)
        V = np.transpose(V)
        R = V.dot(np.transpose(U))

        C_A = np.eye(3)
        C_A = np.vstack(
            [
                np.hstack([C_A, np.transpose(centroid_A) * -1.0]),
                np.array([0, 0, 0, 1]),
            ]
        )

        R_new = np.vstack(
            [np.hstack([R, np.array([[0], [0], [0]])]), np.array([0, 0, 0, 1])]
        )

        C_B = np.eye(3)
        C_B = np.vstack(
            [np.hstack([C_B, np.transpose(centroid_B)]), np.array([0, 0, 0, 1])]
        )

        T = C_B.dot(R_new).dot(C_A)

        super(RigidTransform, self).__init__(T)


def unit_vector(v):
    """Return the unit vector of the vector."""
    return v / norm(v)


def angle(u, v, w=None):
    """Return the angle in radians between two vectors."""
    if w is not None:
        u = u - v
        v = w - v
    c = np.dot(u, v) / norm(u) / norm(v)
    return np.arccos(np.clip(c, -1, 1))


def _create_equivalence_transform(equiv):
    """Compute an equivalence transformation.

    Transforms this compound to another compound's coordinate system.

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
            raise ValueError("Equivalence pair not a 2-tuple")
        if not (
            isinstance(pair[0], Compound) and isinstance(pair[1], Compound)
        ):
            raise ValueError(
                f"Equivalence pair type mismatch: pair[0] is a {pair[0]} "
                f"and pair[1] is a {pair[1]}"
            )

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


def equivalence_transform(
    compound, from_positions, to_positions, add_bond=True
):
    """Compute an affine transformation.

    Maps the from_positions to the respective to_positions, and applies this
    transformation to the compound.

    Parameters
    ----------
    compound : mb.Compound
        The Compound to be transformed.
    from_positions : np.ndarray, shape=(n, 3), dtype=float
        Original positions.
    to_positions : np.ndarray, shape=(n, 3), dtype=float
        New positions.
    """
    warn(
        "The `equivalence_transform` function is being phased out in favor of"
        " `force_overlap`.",
        DeprecationWarning,
    )
    from mbuild.port import Port

    T = None
    if isinstance(from_positions, (list, tuple)) and isinstance(
        to_positions, (list, tuple)
    ):
        equivalence_pairs = zip(from_positions, to_positions)
    elif isinstance(from_positions, Port) and isinstance(to_positions, Port):
        equivalence_pairs, T = _choose_correct_port(
            from_positions, to_positions
        )
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
                warn("Attempting to form bond from port that has no anchor")
            else:
                from_positions.anchor.parent.add_bond(
                    (from_positions.anchor, to_positions.anchor)
                )
                to_positions.anchor.parent.add_bond(
                    (from_positions.anchor, to_positions.anchor)
                )


def _choose_correct_port(from_port, to_port):
    """Chooses the direction when using an equivalence transform on two Ports.

    Each Port object actually contains 2 sets of 4 atoms, either of which can be
    used to make a connection with an equivalence transform. This function
    chooses the set of 4 atoms that makes the anchor atoms not overlap which is
    the intended behavior for most use-cases.

    Parameters
    ----------
    from_port : mb.Port
    to_port : mb.Port

    Returns
    -------
    equivalence_pairs : tuple of Ports, shape=(2,)
        Technically, a tuple of the Ports' sub-Compounds ('up' or 'down') that
        are used to make the correct connection between components.
    """
    # First we try matching the two 'up' ports.
    T1 = _create_equivalence_transform([(from_port["up"], to_port["up"])])
    new_position = T1.apply_to(np.array(from_port.anchor.pos, ndmin=2))

    dist_between_anchors_up_up = norm(new_position[0] - to_port.anchor.pos)

    # Then matching a 'down' with an 'up' port.
    T2 = _create_equivalence_transform([(from_port["down"], to_port["up"])])
    new_position = T2.apply_to(np.array(from_port.anchor.pos, ndmin=2))

    # Determine which transform places the anchors further away from each other.
    dist_between_anchors_down_up = norm(new_position[0] - to_port.anchor.pos)
    difference_between_distances = (
        dist_between_anchors_down_up - dist_between_anchors_up_up
    )

    if difference_between_distances > 0:
        correct_port = from_port["down"]
        T = T2
    else:
        correct_port = from_port["up"]
        T = T1
    return [(correct_port, to_port["up"])], T


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
        raise ValueError("Cannot rotate around a zero vector")
    return Rotation(theta, around).apply_to(coordinates)


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
        raise ValueError("Cannot spin around a zero vector")
    center_pos = np.mean(coordinates, axis=0)
    coordinates -= center_pos
    coordinates = _rotate(coordinates, theta, around)
    coordinates += center_pos
    return coordinates


def x_axis_transform(
    compound, new_origin=None, point_on_x_axis=None, point_on_xy_plane=None
):
    """Move a compound such that the x-axis lies on specified points.

    Parameters
    ----------
    compound : mb.Compound
        The compound to move.
    new_origin : mb.Compound or list-like of size 3,  default=[0.0, 0.0, 0.0]
        Where to place the new origin of the coordinate system.
    point_on_x_axis : mb.Compound or list-like of size 3, default=[1, 0, 0]
        A point on the new x-axis.
    point_on_xy_plane : mb.Compound or list-like of size 3, default=[1, 0, 0]
        A point on the new xy-plane.
    """
    import mbuild as mb

    if new_origin is None:
        new_origin = np.array([0, 0, 0])
    elif isinstance(new_origin, mb.Compound):
        new_origin = new_origin.pos
    elif isinstance(new_origin, (tuple, list, np.ndarray)):
        new_origin = np.asarray(new_origin)
    else:
        raise TypeError(
            "x_axis_transform, y_axis_transform, and z_axis_transform only "
            "accept mb.Compounds, list-like of length 3 or None for the "
            f"new_origin parameter. User passed type: {type(new_origin)}."
        )
    if point_on_x_axis is None:
        point_on_x_axis = np.array([1.0, 0.0, 0.0])
    elif isinstance(point_on_x_axis, mb.Compound):
        point_on_x_axis = point_on_x_axis.pos
    elif isinstance(point_on_x_axis, (list, tuple, np.ndarray)):
        point_on_x_axis = np.asarray(point_on_x_axis)
    else:
        raise TypeError(
            "x_axis_transform, y_axis_transform, and z_axis_transform only "
            "accept mb.Compounds, list-like of size 3, or None for the "
            "point_on_x_axis parameter. User passed type: "
            f"{point_on_x_axis}."
        )
    if point_on_xy_plane is None:
        point_on_xy_plane = np.array([1.0, 1.0, 0.0])
    elif isinstance(point_on_xy_plane, mb.Compound):
        point_on_xy_plane = point_on_xy_plane.pos
    elif isinstance(point_on_xy_plane, (list, tuple, np.ndarray)):
        point_on_xy_plane = np.asarray(point_on_xy_plane)
    else:
        raise TypeError(
            "x_axis_transform, y_axis_transform, and z_axis_transform only "
            "accept mb.Compounds, list-like of size 3, or None for the "
            "point_on_xy_plane parameter. User passed type: "
            f"{point_on_xy_plane}."
        )

    atom_positions = compound.xyz_with_ports
    transform = AxisTransform(
        new_origin=new_origin,
        point_on_x_axis=point_on_x_axis,
        point_on_xy_plane=point_on_xy_plane,
    )
    atom_positions = transform.apply_to(atom_positions)
    compound.xyz_with_ports = atom_positions


def y_axis_transform(
    compound, new_origin=None, point_on_y_axis=None, point_on_xy_plane=None
):
    """Move a compound such that the y-axis lies on specified points.

    Parameters
    ----------
    compound : mb.Compound
        The compound to move.
    new_origin : mb.Compound or like-like of size 3, default=[0, 0, 0]
        Where to place the new origin of the coordinate system.
    point_on_y_axis : mb.Compound or list-like of size 3, default=[0, 1, 0]
        A point on the new y-axis.
    point_on_xy_plane : mb.Compound or list-like of size 3, default=[0, 1, 0]
        A point on the new xy-plane.
    """
    x_axis_transform(
        compound,
        new_origin=new_origin,
        point_on_x_axis=point_on_y_axis,
        point_on_xy_plane=point_on_xy_plane,
    )
    compound.rotate(theta=np.pi / 2, around=(0, 0, 1))


def z_axis_transform(
    compound, new_origin=None, point_on_z_axis=None, point_on_zx_plane=None
):
    """Move a compound such that the z-axis lies on specified points.

    Parameters
    ----------
    compound : mb.Compound
        The compound to move.
    new_origin : mb.Compound or list-like of size 3, default=[0, 0, 0]
        Where to place the new origin of the coordinate system.
    point_on_z_axis : mb.Compound or list-like of size 3, default=[0, 0, 1]
        A point on the new z-axis.
    point_on_zx_plane : mb.Compound or list-like of size 3, default=[0, 0, 1]
        A point on the new xz-plane.
    """
    x_axis_transform(
        compound,
        new_origin=new_origin,
        point_on_x_axis=point_on_z_axis,
        point_on_xy_plane=point_on_zx_plane,
    )
    compound.rotate(theta=np.pi * 1.5, around=(0, 1, 0))
