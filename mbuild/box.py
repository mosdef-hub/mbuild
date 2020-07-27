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

    @classmethod
    def from_lengths_angles(cls, lengths, angles):
        pass

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
        pass

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

class _BoxArray(np.ndarray):
    """Subclass of np.ndarry specifically for mb.Box

    This subclass is meant to be used internally to store Box attribute array.
    This subclass is modified so that its __setitem__ method is rerouted to the
    corresponding setter method.

    Parameters
    ----------
    array : array-like object
        This can be tuple, list, or any array-like object that can be usually
        passed to np.array
    var : str
        Corresponding Box's attributes like "maxs", "mins", "lengths", "angles"
    box : mb.Box
        This is the Box contains this attribute (one level up of this array)
    """
    def __new__(cls, array, var=None, box=None, dtype=np.float):
        _array = np.asarray(array, dtype).view(cls)
        _array.var = var
        _array.box = box
        return _array

    def __setitem__(self, key, val):
        array = list(self)
        array[key] = val
        if self.var == "maxs":
            msg = "Given max value is less than box's min value"
            assert val >= self.box.mins[key], msg
            self.box.maxs = array
        elif self.var == "mins":
            msg = "Given min value is more than box's max value"
            assert val <= self.box.maxs[key], msg
            self.box.mins = array
        elif self.var == "lengths":
            msg = "Given length value is negative"
            assert val >= 0, msg
            self.box.lengths = array
        else:
            self.box.angles = array
