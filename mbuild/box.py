from warnings import warn

import numpy as np


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
            self._mins = BoxArray(array=mins, var="mins", box=self)
            self._maxs = BoxArray(array=maxs, var="maxs", box=self)
            self._lengths = BoxArray(array=(self.maxs - self.mins), var="lengths", box=self)
        else:
            if mins is not None or maxs is not None:
                warn(
                    "Provided `lengths` and `mins` and/or `maxs`. Only `lengths` "
                    "is being used. You provided: "
                    "lengths={} mins={} maxs={}".format(lengths, mins, maxs)
                )
            self._mins = BoxArray(array=(0,0,0), var="mins", box=self)
            self._maxs = BoxArray(array=lengths, var="maxs", box=self)
            self._lengths = BoxArray(array=lengths, var="lenghts", box=self)
        if angles is None:
            angles = BoxArray(array=(90.0, 90.0, 90.0), var="angles", box=self)
        elif isinstance(angles, (list, np.ndarray)):
            angles = BoxArray(array=angles, var="anglese", box=self)
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
        self._mins = BoxArray(array=mins, var="mins", box=self)
        self._lengths = BoxArray(array=(self.maxs - self.mins), var="lengths", box=self)

    @maxs.setter
    def maxs(self, maxs):
        maxs = np.array(maxs, dtype=np.float)
        assert maxs.shape == (3, )
        self._maxs = BoxArray(array=maxs, var="maxs", box=self)
        self._lengths = BoxArray(array= (self.maxs - self.mins), var="lengths", box=self)

    @lengths.setter
    def lengths(self, lengths):
        lengths = np.array(lengths, dtype=np.float)
        assert lengths.shape == (3, )
        self._maxs = BoxArray(array=(self.maxs + (0.5*lengths - 0.5*self.lengths)), var="maxs", box=self)
        self._mins = BoxArray(array=(self.mins - (0.5*lengths - 0.5*self.lengths)), var="mins", box=self)
        self._lengths = BoxArray(array=lengths, var="lengths", box=self, dtype=np.float)

    @angles.setter
    def angles(self, angles):
        angles = np.array(angles, dtype=np.float)
        assert angles.shape == (3, )
        self._angles = angles

    def __repr__(self):
        return "Box(mins={}, maxs={}, angles={})".format(self.mins, self.maxs, self.angles)

class BoxArray(np.ndarray):
    """Subclass of np.ndarry specifically for mb.Box

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
            self.box.maxs = array
        elif self.var == "mins":
            self.box.mins = array
        elif self.var == "lengths":
            self.box.lengths = array
        else:
            self.angles = array
