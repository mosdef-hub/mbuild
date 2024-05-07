"""mBuild box module."""

from warnings import warn

import numpy as np

from mbuild.exceptions import MBuildError

__all__ = ["Box"]


class Box(object):
    """A box representing the bounds of the system.

    Parameters
    ----------
    lengths : list-like, shape=(3,), dtype=float
        Lengths of the edges of the box.
    angles : list-like, shape=(3,), dtype=float, default=None
        Angles (in degrees) that define the tilt of the edges of the box. If
        None is given, angles are assumed to be [90.0, 90.0, 90.0]. These are
        also known as alpha, beta, gamma in the crystallography community.
    precision : int, optional, default=None
        Control the precision of the floating point representation of box
        attributes. If none provided, the default is 6 decimals.

    Attributes
    ----------
    vectors : np.ndarray, shape=(3,3), dtype=float
        Vectors that define the parallelepiped (Box).
    lengths : tuple, shape=(3,), dtype=float
        Lengths of the box in x,y,z
    angles : tuple, shape=(3,), dtype=float
        Angles defining the tilt of the box.
    Lx : float
        Length of the Box in the x dimension
    Ly : float
        Length of the Box in the y dimension
    Lz : float
        Length of the Box in the z dimension
    xy : float
        Tilt factor needed to displace an orthogonal box's xy face to its
        parallelepiped structure.
    xz : float
        Tilt factor needed to displace an orthogonal box's xz face to its
        parallelepiped structure.
    yz : float
        Tilt factor needed to displace an orthogonal box's yz face to its
        parallelepiped structure.
    precision : int
        Precision of the floating point numbers when accessing values.

    Notes
    -----
    Box vectors are expected to be provided in row-major format.
    """

    def __init__(self, lengths, angles=None, precision=None):
        if precision is not None:
            self._precision = int(precision)
        else:
            self._precision = 6

        if angles is None:
            angles = [90.0, 90.0, 90.0]

        self._vectors = _lengths_angles_to_vectors(
            lengths=lengths, angles=angles, precision=self.precision
        )
        (Lx, Ly, Lz, xy, xz, yz) = self._from_vecs_to_lengths_tilt_factors()
        self._Lx = Lx
        self._Ly = Ly
        self._Lz = Lz
        self._xy = xy
        self._xz = xz
        self._yz = yz

    @classmethod
    def from_lengths_angles(cls, lengths, angles, precision=None):
        """Generate a box from lengths and angles."""
        return cls(lengths=lengths, angles=angles, precision=precision)

    @classmethod
    def from_uvec_lengths(cls, uvec, lengths, precision=None):
        """Generate a box from unit vectors and lengths."""
        uvec = np.asarray(uvec)
        uvec.reshape(3, 3)

        if not np.allclose(np.linalg.norm(uvec, axis=1), 1.0):
            raise MBuildError(
                "Unit vector magnitudes provided are not close to 1.0, "
                f"magnitudes: {np.linalg.norm(uvec, axis=1)}"
            )

        lengths = np.asarray(lengths)
        lengths.reshape(1, 3)
        _validate_box_vectors(uvec)
        scaled_vec = (uvec.T * lengths).T
        (alpha, beta, gamma) = _calc_angles(scaled_vec)

        return cls(
            lengths=lengths, angles=(alpha, beta, gamma), precision=precision
        )

    @classmethod
    def from_mins_maxs_angles(cls, mins, maxs, angles, precision=None):
        """Generate a box from min/max distance calculations and angles."""
        (x_min, y_min, z_min) = mins
        (x_max, y_max, z_max) = maxs
        lengths = (x_max - x_min, y_max - y_min, z_max - z_min)
        return cls(lengths=lengths, angles=angles, precision=precision)

    @classmethod
    def from_vectors(cls, vectors, precision=None):
        """Generate a box from box vectors."""
        vectors = _validate_box_vectors(vectors)
        (alpha, beta, gamma) = _calc_angles(vectors)
        v1 = vectors[0, :]
        v2 = vectors[1, :]
        v3 = vectors[2, :]

        Lx = np.linalg.norm(v1)
        Ly = np.linalg.norm(v2)
        Lz = np.linalg.norm(v3)
        lengths = (Lx, Ly, Lz)
        return cls(
            lengths=lengths, angles=(alpha, beta, gamma), precision=precision
        )

    @classmethod
    def from_lengths_tilt_factors(
        cls, lengths, tilt_factors=None, precision=None
    ):
        """Generate a box from box lengths and tilt factors."""
        (Lx, Ly, Lz) = lengths
        if tilt_factors is None:
            (xy, xz, yz) = (0.0, 0.0, 0.0)
        else:
            (xy, xz, yz) = tilt_factors

        vecs = np.asarray(
            [[Lx, 0.0, 0.0], [Ly * xy, Ly, 0.0], [Lz * xz, Lz * yz, Lz]]
        )
        (alpha, beta, gamma) = _calc_angles(vecs)
        return cls(
            lengths=lengths, angles=[alpha, beta, gamma], precision=precision
        )

    @classmethod
    def from_lo_hi_tilt_factors(cls, lo, hi, tilt_factors, precision=None):
        """Generate a box from a lo, hi convention and tilt factors."""
        (xlo, ylo, zlo) = lo
        (xhi, yhi, zhi) = hi
        (xy, xz, yz) = tilt_factors

        xlo_bound = xlo + min([0.0, xy, xz, xy + xz])
        xhi_bound = xhi + max([0.0, xy, xz, xy + xz])
        ylo_bound = ylo + min([0.0, yz])
        yhi_bound = yhi + max([0.0, yz])

        lengths = [xhi_bound - xlo_bound, yhi_bound - ylo_bound, zhi - zlo]
        return cls.from_lengths_tilt_factors(
            lengths=lengths, tilt_factors=tilt_factors
        )

    @property
    def vectors(self):
        """Box representation as a 3x3 matrix."""
        return self._vectors

    @property
    def box_parameters(self):
        """Lengths and tilt factors of the box."""
        return self.Lx, self.Ly, self.Lz, self.xy, self.xz, self.xy

    @property
    def Lx(self):
        """Length in the x direction."""
        return round(self._Lx, self.precision)

    @property
    def Ly(self):
        """Length in the y direction."""
        return round(self._Ly, self.precision)

    @property
    def Lz(self):
        """Length in the z direction."""
        return round(self._Lz, self.precision)

    @property
    def lengths(self):
        """Lengths of the box."""
        return self.Lx, self.Ly, self.Lz

    @property
    def xy(self):
        """Tilt factor xy of the box."""
        return round(self._xy, self.precision)

    @property
    def xz(self):
        """Tilt factor xz of the box."""
        return round(self._xz, self.precision)

    @property
    def yz(self):
        """Tilt factor yz of the box."""
        return round(self._yz, self.precision)

    @property
    def tilt_factors(self):
        """Return the 3 tilt_factors (xy, xz, yz) of the box."""
        return self.xy, self.xz, self.yz

    @property
    def angles(self):
        """Angles defining the tilt of the box (alpha, beta, gamma)."""
        (alpha, beta, gamma) = self._get_angles()
        alpha = round(alpha, self.precision)
        beta = round(beta, self.precision)
        gamma = round(gamma, self.precision)
        return alpha, beta, gamma

    @property
    def precision(self):
        """Amount of decimals to represent floating point values."""
        return self._precision

    @precision.setter
    def precision(self, value):
        """Decimal point precision, if None use 16, else cast as int."""
        if not value:
            precision = 16
        else:
            precision = int(value)
        self._precision = precision

    @property
    def bravais_parameters(self):
        """Return the Box representation as Bravais lattice parameters.

        Based on the box vectors, return the parameters to describe the box in
        terms of the Bravais lattice parameters:

            a,b,c = the edges of the Box
            alpha, beta, gamma = angles(tilt) of the parallelepiped, in degrees

        Returns
        -------
        parameters : tuple of floats,
            (a, b, c, alpha, beta, gamma)
        """
        (alpha, beta, gamma) = self.angles
        (Lx, Ly, Lz) = self.lengths
        return Lx, Ly, Lz, alpha, beta, gamma

    def __repr__(self):
        """Return a string representation of the box."""
        (Lx, Ly, Lz, xy, xz, yz) = self.box_parameters
        format_precision = f".{self._precision}f" if self._precision else ""
        desc = (
            f"Box: Lx={Lx:{format_precision}}, "
            f"Ly={Ly:{format_precision}}, "
            f"Lz={Lz:{format_precision}}, "
            f"xy={xy:{format_precision}}, "
            f"xz={xz:{format_precision}}, "
            f"yz={yz:{format_precision}}, "
        )
        return desc

    def _from_vecs_to_lengths_tilt_factors(self):
        # vectors should already be aligned by _normalize_box
        v = np.zeros((3, 3))
        v[0, :] = self._vectors[0, :]
        v[1, :] = self._vectors[1, :]
        v[2, :] = self._vectors[2, :]

        Lx = np.sqrt(np.dot(v[0], v[0]))
        a2x = np.dot(v[0], v[1]) / Lx
        Ly = np.sqrt(np.dot(v[1], v[1]) - a2x * a2x)
        xy = a2x / Ly
        v0xv1 = np.cross(v[0], v[1])
        v0xv1mag = np.sqrt(np.dot(v0xv1, v0xv1))
        Lz = np.dot(v[2], v0xv1) / v0xv1mag
        a3x = np.dot(v[0], v[2]) / Lx
        xz = a3x / Lz
        yz = (np.dot(v[1], v[2]) - a2x * a3x) / (Ly * Lz)

        len_x = np.sqrt(np.dot(v[0], v[0]))
        len_y = np.sqrt(np.dot(v[1], v[1]))
        len_z = np.sqrt(np.dot(v[2], v[2]))
        return len_x, len_y, len_z, xy, xz, yz

    def _get_angles(self):
        return _calc_angles(self.vectors)


def _validate_box_vectors(box_vectors):
    """Determine if the vectors are in the convention we use.

    This method will parse the provided box vectors, determine if the vectors
    follow the conventions the Box class adheres to. In this case:

    1. It is a 3x3 matrix that can be coerced into a numpy array of floats
    2. Vectors are in a right-handed basis (determinant of matrix is +)
    3. The first vector is aligned along the [1,0,0] direction
    4. The second vector in aligned along the xy plane
    5. The third vector can align freely in the x,y, and +z direction

    If the vectors are not right-handed, a warning will be raised, and the
    vectors will be transformed into the right-handed coordinate system.

    If the three vectors are not following conventions 3-5, the matrix will be
    transformed to comply with them, and will also raise a warning.
    """
    vecs = np.asarray(box_vectors, dtype=np.float64)
    vecs.reshape(3, 3)

    return _normalize_box(vecs)


def _lengths_angles_to_vectors(lengths, angles, precision):
    (a, b, c) = lengths
    (alpha, beta, gamma) = np.deg2rad(angles)
    cos_a = np.clip(np.cos(alpha), -1.0, 1.0)
    cos_b = np.clip(np.cos(beta), -1.0, 1.0)
    cos_g = np.clip(np.cos(gamma), -1.0, 1.0)

    sin_a = np.clip(np.sin(alpha), -1.0, 1.0)
    sin_b = np.clip(np.sin(beta), -1.0, 1.0)
    sin_g = np.clip(np.sin(gamma), -1.0, 1.0)
    a_vec = np.asarray([a, 0.0, 0.0])

    b_x = b * cos_g
    b_y = b * sin_g
    b_vec = np.asarray([b_x, b_y, 0.0])

    c_x = c * cos_b
    c_cos_y_term = (cos_a - (cos_b * cos_g)) / sin_g
    c_y = c * c_cos_y_term
    c_z = c * np.sqrt(1 - np.square(cos_b) - np.square(c_cos_y_term))
    c_vec = np.asarray([c_x, c_y, c_z])
    box_vectors = np.asarray((a_vec, b_vec, c_vec))
    box_vectors.reshape(3, 3)
    _validate_box_vectors(box_vectors=box_vectors)
    return box_vectors.round(precision)


def _normalize_box(vectors):
    """Align the box matrix into a right-handed coordinate frame.

    NOTE: This assumes that the matrix is in a row-major format.

    NOTE: Inspiration and logic are from the Glotzer group package, Garnett;
    which is provided under a BSD 3-clause License.
    For additional information, refer to the License file provided with this
    package.
    """
    det = np.linalg.det(vectors)
    if np.isclose(det, 0.0, atol=1e-5):
        raise MBuildError(
            "The vectors to define the box are co-linear, this does not form a "
            f"3D region in space.\n Box vectors evaluated: {vectors}"
        )
    if det < 0.0:
        warn(
            "Box vectors provided for a left-handed basis, these will be "
            "transformed into a right-handed basis automatically."
        )

    # transpose to column-major for the time being
    Q, R = np.linalg.qr(vectors.T)

    # left or right handed: det<0 left, >0, right
    sign = np.linalg.det(Q)
    R = R * sign

    signs = np.diag(
        np.diag(np.where(R < 0, -np.ones(R.shape), np.ones(R.shape)))
    )
    transformed_vecs = R.dot(signs)
    return _reduced_form_vectors(transformed_vecs.T)


def _reduced_form_vectors(box_vectors):
    """Get reduced vectors from vectors.

    Adapted from HOOMD-Blue's documentation on periodic boundary conditions:
    https://hoomd-blue.readthedocs.io/en/stable/box.html
    """
    v1 = box_vectors[0, :]
    v2 = box_vectors[1, :]
    v3 = box_vectors[2, :]

    lx = np.linalg.norm(v1)
    a_2x = np.dot(v1, v2) / lx
    ly = np.sqrt(np.dot(v2, v2) - a_2x * a_2x)
    xy = a_2x / ly
    v1_x_v2 = np.cross(v1, v2)
    lz = np.dot(v3, (v1_x_v2 / np.linalg.norm(v1_x_v2)))
    a_3x = np.dot(v1, v3) / lx
    xz = a_3x / lz
    yz = (np.dot(v2, v3) - a_2x * a_3x) / (ly * lz)

    reduced_vecs = np.asarray(
        [[lx, 0.0, 0.0], [xy * ly, ly, 0.0], [xz * lz, yz * lz, lz]]
    )

    return reduced_vecs


def _calc_angles(vectors):
    """Calculate the angles between the vectors that define the box.

    Calculates the angles alpha, beta, and gamma from the Box object
    attribute box_vectors, rounded to 'precision' number of decimal points.
    """
    vector_magnitudes = np.linalg.norm(vectors, axis=1)
    v = np.zeros((3, 3))
    v[0, :] = vectors[0, :]
    v[1, :] = vectors[1, :]
    v[2, :] = vectors[2, :]

    a_dot_b = np.dot(vectors[0], vectors[1])
    b_dot_c = np.dot(vectors[1], vectors[2])
    a_dot_c = np.dot(vectors[0], vectors[2])

    alpha_raw = b_dot_c / (vector_magnitudes[1] * vector_magnitudes[2])
    beta_raw = a_dot_c / (vector_magnitudes[0] * vector_magnitudes[2])
    gamma_raw = a_dot_b / (vector_magnitudes[0] * vector_magnitudes[1])

    (alpha, beta, gamma) = np.rad2deg(
        np.arccos(np.clip([alpha_raw, beta_raw, gamma_raw], -1.0, 1.0))
    )

    return alpha, beta, gamma
