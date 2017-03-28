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
    def __init__(self, lengths=None, mins=None, maxs=None):
        if lengths is not None:
            assert mins is None and maxs is None
            self._mins = np.array([0.0, 0.0, 0.0])
            self._maxs = np.array(lengths)
            self._lengths = np.array(lengths)
        elif maxs is not None:
            assert mins is not None and lengths is None
            self._mins = np.array(mins)
            self._maxs = np.array(maxs)
            self._lengths = self.maxs - self.mins
        else:
            raise ValueError("Either provide `lengths` or `mins` and `maxs`."
                             "You provided: lengths={} mins={} maxs={}".format(lengths, mins, maxs))

    @property
    def mins(self):
        return self._mins

    @property
    def maxs(self):
        return self._maxs

    @property
    def lengths(self):
        return self._lengths

    @mins.setter
    def mins(self, mins):
        assert mins.shape == (3, )
        self._mins = mins
        self._lengths = self.maxs - self.mins

    @maxs.setter
    def maxs(self, maxes):
        assert maxes.shape == (3, )
        self._maxs = maxes
        self._lengths = self.maxs - self.mins

    @lengths.setter
    def lengths(self, lengths):
        assert lengths.shape == (3, )
        self._maxs += 0.5*lengths - 0.5*self.lengths
        self._mins -= 0.5*lengths - 0.5*self.lengths
        self._lengths = lengths

    def __repr__(self):
        return "Box(mins={}, maxs={})".format(self.mins, self.maxs)
