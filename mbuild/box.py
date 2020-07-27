from warnings import warn

import numpy as np

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

        self.box_vectors = box_vectors

    @classmethod
    def from_lengths_angles(cls, lengths, angles):
        box_vectors = _lengths_angles_to_vectors(lengths, angles)
        _validate_box_vectors(box_vectors)

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
        pass

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

    @mins.setter
    def mins(self, mins):
        mins = np.array(mins, dtype=np.float)
        assert mins.shape == (3, )
        assert all(mins <= self.maxs), "Given mins is greater than maxs"
        self._mins = _BoxArray(array=mins, var="mins", box=self)
        self._lengths = _BoxArray(array=(self.maxs - self.mins), var="lengths", box=self)

    @maxs.setter
    def maxs(self, maxs):
        maxs = np.array(maxs, dtype=np.float)
        assert maxs.shape == (3, )
        assert all(maxs >= self.mins), "Given maxs is less than mins"
        self._maxs = _BoxArray(array=maxs, var="maxs", box=self)
        self._lengths = _BoxArray(array= (self.maxs - self.mins), var="lengths", box=self)

    @lengths.setter
    def lengths(self, lengths):
        if isinstance(lengths, int) or isinstance(lengths, float):
            lengths = np.array(lengths*np.ones(3), dtype=np.float)
        else:
            lengths = np.array(lengths, dtype=np.float)
        assert lengths.shape == (3, )
        assert all(lengths >= 0), "Given lengths are negative"
        self._maxs = _BoxArray(array=(self.maxs + (0.5*lengths - 0.5*self.lengths)), var="maxs", box=self)
        self._mins = _BoxArray(array=(self.mins - (0.5*lengths - 0.5*self.lengths)), var="mins", box=self)
        self._lengths = _BoxArray(array=lengths, var="lengths", box=self, dtype=np.float)

    @angles.setter
    def angles(self, angles):
        angles = np.array(angles, dtype=np.float)
        assert angles.shape == (3, )
        self._angles = _BoxArray(array=angles, var="angles", box=self, dtype=np.float)

    def __repr__(self):
        return "Box(mins={}, maxs={}, angles={})".format(self.mins, self.maxs, self.angles)


def _validate_box_vectors(box_vectors):
    pass

def _lengths_angles_to_vectors(lengths, angles):
    (a, b, c) = lengths
    (alpha, beta, gamma) = np.deg2rad(angles)

    a_vec = np.asarray([a,0.0,0.0],)

    b_x = b * np.cos(gamma)
    b_y = b * np.sin(gamma)
    b_vec = np.asarray([b_x, b_y, 0.0],))

    c_x = c * np.cos(beta)
    c_cos_y_term = ((np.cos(alpha) - (np.cos(beta) * np.cos(gamma))) / np.sin(gamma))
    c_y = c * c_cos_y_term
    c_z = c * np.sqrt(1 - np.square(np.cos(beta)) - np.square(c_cos_y_term))
    c_vec = np.asarray([c_x, c_y, c_z])

    box_vectors = np.asarray([a_vec],[b_vec],[c_vec])
    box_vectors.reshape(3,3)
    return box_vectors
