from warnings import warn

import numpy as np

__all__ = ['Box']

class Box(object):
    """A box representing the bounds of the system.

    Attributes
    ----------
    mins : np.ndarray, shape=(3,), dtype=float
        Minimum x, y, z coordinates.
    maxs : np.ndarray, shape=(3,), dtype=float
        Maximum x, y, z coordinates.
    lengths : np.ndarray, shape(3,), dtype=float
        Box length in x, y and z directions.

    """
    def __init__(self, lengths=None, mins=None, maxs=None, angles=None):
        if lengths is None:
            if mins is None or maxs is None:
                raise ValueError(
                    "Either provide `lengths` or `mins` and `maxs`. "
                    "You provided: "
                    "lengths={} mins={} maxs={}".format(lengths, mins, maxs)
                )
            mins = np.array(mins, dtype=np.float)
            maxs = np.array(maxs, dtype=np.float)
            assert mins.shape == (3, ), "Given mins have wrong dimensions"
            assert maxs.shape == (3, ), "Given maxs have wrong dimensions"
            assert all(mins <= maxs), "Given mins are greater than maxs"
            self._mins = _BoxArray(array=mins, var="mins", box=self)
            self._maxs = _BoxArray(array=maxs, var="maxs", box=self)
            self._lengths = _BoxArray(array=(self.maxs - self.mins), var="lengths", box=self)
        else:
            if mins is not None or maxs is not None:
                warn(
                    "Provided `lengths` and `mins` and/or `maxs`. Only `lengths` "
                    "is being used. You provided: "
                    "lengths={} mins={} maxs={}".format(lengths, mins, maxs)
                )
            if isinstance(lengths, int) or isinstance(lengths, float):
                lengths = np.array(lengths*np.ones(3), dtype=np.float)
            else:
                lengths = np.array(lengths, dtype=np.float)
            assert lengths.shape == (3, )
            assert all(lengths >= 0), "Given lengths are negative"
            self._mins = _BoxArray(array=(0,0,0), var="mins", box=self)
            self._maxs = _BoxArray(array=lengths, var="maxs", box=self)
            self._lengths = _BoxArray(array=lengths, var="lengths", box=self)
        if angles is None:
            angles = _BoxArray(array=(90.0, 90.0, 90.0), var="angles", box=self)
        elif isinstance(angles, (list, np.ndarray)):
            angles = _BoxArray(array=angles, var="angles", box=self)
        self._angles = angles

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
    This subclass is modified so that its __setitem__ method is reroute to the
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
