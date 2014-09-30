__author__ = 'sallai'

import numpy as np


class Box(object):
    """A box representing the bounds of the system.

    Attributes:
       mins (np.ndarray of floats): Minimum x, y, z coordinates.
       maxs (np.ndarray of floats): Maximum x, y, z coordinates.
       lengths (np.ndarray of floats): Box length in x, y and z directions. 

    """

    def __init__(self, lengths=None, mins=None, maxs=None):
        """Construct a simulation box in one of two ways.

        1) Specifying the box length sets the mins to 0.0 and the maxs to
        lengths.
        2) Specifying the maxs and mins sets the lengths to maxs-mins.

        Adjusting the length of a box adjusts `maxs` to reflect the change.

        Args:
            mins (np.ndarray of floats, optional): Minimum x, y, z
                coordinates.
            maxs (np.ndarray of floats, optional): Maximum x, y, z coordinates.
            lengths (np.ndarray of floats, optional): Box length in x, y and z
                directions. 

        """
        if lengths is not None:
            assert(mins is None and maxs is None)
            self._mins = np.array([0.0, 0.0, 0.0])
            self._maxs = np.array(lengths)
            self._lengths = np.array(lengths)

        elif maxs is not None:
            assert(mins is not None and lengths is None)
            self._mins = np.array(mins)
            self._maxs = np.array(maxs)
            self._lengths = self.maxs - self.mins

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
