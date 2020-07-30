from warnings import warn

import numpy as np

from mbuild.coordinate_transform import AxisTransform, ChangeOfBasis
from mbuild.exceptions import MBuildError

__all__ = ['Box']

class Box(object):
    """A box representing the bounds of the system.

    Parameters
    ----------
    box_vectors : np.ndarray, shape=(3,3), dtype=float
        Vectors that define a right-handed parallelpiped (Box).


    Attributes
    ----------
    box_vectors : np.ndarray, shape=(3,3), dtype=float
        Vectors that define the parallelpiped (Box).
    Lx, Ly, Lz : float
        Lengths of the Box in the x,y,z dimensions
    xy,xz,yz : float
        Tilt factors needed to displace an orthogonal box to its parallelpiped structure.
    """
    def __init__(self, box_vectors=None):
        try:
            _validate_box_vectors(box_vectors)
        except:
            pass

        self._box_vectors = box_vectors

    @classmethod
    def from_lengths_angles(cls, lengths, angles):
        box_vectors = _lengths_angles_to_vectors(lengths, angles)

        return Box(box_vectors=box_vectors)

    @classmethod
    def from_uvec_lengths(cls, uvec, lengths):
        uvec = np.asarray(uvec)
        uvec.reshape(3,3)
        assert uvec.shape == (3,3), f"Expected a 3x3 matrix, was provided {uvec.shape}."

        lengths = np.asarray(lengths)
        lengths.reshape(1,3)
        _validate_box_vectors(uvec)
        scaled_vec = (uvec.T * lengths).T

        return Box(box_vectors=scaled_vec)

    @classmethod
    def from_mins_maxs_angles(cls, mins, maxs, angles):
        (x_min, y_min, z_min) = mins
        (x_max, y_max, z_max) = maxs
        lengths = (x_max-x_min, y_max-y_min, z_max-z_min)
        box_vectors = _lengths_angles_to_vectors(lengths, angles)
        return Box(box_vectors=box_vectors)

    @classmethod
    def from_lengths_tilt_factors(cls, lengths, tilt_factors):
        pass

    @classmethod
    def from_lo_hi_tilt_factors(cls, lo, hi, tilt_factors):
        (xlo,ylo,zlo) = lo
        (xhi,yhi,zhi) = hi
        (xy,xz,yz) = tilt_factors

        v1 = np.asarray([xhi - xlo, 0.0, 0.0])
        v2 = np.asarray([xy, yhi - ylo, 0.0])
        v3 = np.asarray([xz, yz, zhi - zlo])
        box_vectors = np.asarray([v1],[v2],[v3])
        box_vectors.reshape(3,3)
        _validate_box_vectors(box_vectors)

        return Box(box_vectors=box_vectors)

    @classmethod
    def from_lengths_mins_angles(cls, lengths, mins, angles):
        pass

    @classmethod
    def from_lengths_maxs_angles(cls, lengths, maxs, angles):
        pass

    @property
    def box_vectors(self):
        return self._box_vectors

    #NOTE: we might not want setters at all, just make a new box?
    @box_vectors.setter
    def box_vectors(self, box_vectors):
        self._box_vectors = box_vectors

    @property
    def mins(self):
        return self._mins

    @property
    def maxs(self):
        return self._maxs

    @property
    def lengths(self):
        return self._lengths

    @property
    def angles(self):
        return self._angles

    def __repr__(self):
        return "Box(mins={}, maxs={}, angles={})".format(self.mins, self.maxs, self.angles)


def _validate_box_vectors(box_vectors):
    """Determine if the vectors are in the convention we use.

    This method will parse the inputted box vectors, determine if the
    vectors follow the conventions the Box class adheres to. In this case:

    1. It is a 3x3 matrix that can be coerced into a numpy array of floats
    2. Vectors are in a right-handed basis (determinant of matrix is +)
    3. The first vector is aligned along the [1,0,0] direction
    4. The second vector in aligned along the xy plane
    5. The third vector can align freely in the x,y, and +z direction

    If the vectors are not right-handed, a warning will be raised,
    and the vectors will be transformed into the right-handed coordinate
    system.

    If the three vectors are not following conventions 3-5, the matrix will
    be transformed to comply with them, and will also raise a warning.
    """
    vecs = np.asarray(box_vectors, dtype=np.float)
    if vecs.shape == (3,3):
        pass
    else:
        vecs.reshape(3,3)

    #verify handedness, change handedness if needed
    basis_rh = np.asarray([[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]])
    origin = np.asarray([0.0, 0.0, 0.0])
    vecs = _validate_handedness(basis_rh, origin, vecs)

    # align axes to standard
    vecs = _transform_axis_to_standard(vectors=vecs,
                                       origin=origin,
                                       point_on_x_axis=np.array([1.0, 0.0, 0.0]),
                                       point_on_xy_plane=np.array([1.0, 1.0, 0.0]))
    return vecs

def _lengths_angles_to_vectors(lengths, angles):
    (a, b, c) = lengths
    (alpha, beta, gamma) = np.deg2rad(angles)

    a_vec = np.asarray([a,0.0,0.0],)

    b_x = b * np.cos(gamma)
    b_y = b * np.sin(gamma)
    b_vec = np.asarray([b_x, b_y, 0.0],)

    c_x = c * np.cos(beta)
    c_cos_y_term = ((np.cos(alpha) - (np.cos(beta) * np.cos(gamma))) / np.sin(gamma))
    c_y = c * c_cos_y_term
    c_z = c * np.sqrt(1 - np.square(np.cos(beta)) - np.square(c_cos_y_term))
    c_vec = np.asarray([c_x, c_y, c_z])

    box_vectors = np.asarray([a_vec],[b_vec],[c_vec])
    box_vectors.reshape(3,3)
    _validate_box_vectors(box_vectors=box_vectors)
    return box_vectors

def _validate_handedness(basis, origin, vectors):
    """If the vectors are already right handed, do nothing, else, change basis.
    """
    det = np.linalg.det(vectors)
    if det > 0:
        return vectors
    elif det < 0:
        # left handed basis, convert to right-handed
        warn("Box vectors provided for a left-handed basis, these will "
             "be transformed into a right-handed basis automatically.")
        basis_rh = np.asarray([[1, 0, 0],
                               [0, 1, 0],
                               [0, 0, 1]])
        basis_change = ChangeOfBasis(basis=basis_rh, origin=origin)
        vecs = basis_change.apply_to(vectors)
        assert np.linalg.det(vecs) > 0
        return vecs
    else:
        raise MBuildError("The vectors to define the box are co-linear, "
                          "this does not form a 3D region in space.\n"
                          f"Box vectors evaluated: {vecs}")


def _transform_axis_to_standard(vectors, origin, point_on_x_axis, point_on_xy_plane):
    transform_axis = AxisTransform(new_origin=origin, point_on_x_axis=point_on_x_axis,
                                   point_on_xy_plane=point_on_xy_plane)
    return transform_axis.apply_to(vectors)
