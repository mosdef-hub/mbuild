__author__ = 'sallai'

import numpy as np

class Box(object):
    """A box representing the bounds of the system.

    Attributes:
       mins (np.ndarray of floats): Minimum x, y, z coordinates.
       maxes (np.ndarray of floats): Maximum x, y, z coordinates. 
       lengths (np.ndarray of floats): Box length in x, y and z directions. 

    """

    def __init__(self, lengths=None, mins=None, maxes=None):
        """Construct a simulation box in one of two ways.

        1) Specifying the box length sets the mins to 0.0 and the maxes to
        lengths.
        2) Specifying the maxes and mins sets the lengths to maxes-mins.

        Adjusting the length of a box adjusts `maxes` to reflect the change.

        Args:
            mins (np.ndarray of floats, optional): Minimum x, y, z
                coordinates.
            maxes (np.ndarray of floats, optional): Maximum x, y, z coordinates. 
            lengths (np.ndarray of floats, optional): Box length in x, y and z
                directions. 

        """
        if lengths is not None:
            assert(mins is None and maxes is None)
            self._mins = np.array([0.0, 0.0, 0.0])
            self._maxes = np.array(lengths)
            self._lengths = np.array(lengths)

        elif maxes is not None:
            assert(mins is not None and lengths is None)
            self._mins = np.array(mins)
            self._maxes = np.array(maxes)
            self._lengths = self.maxes - self.mins

    @property
    def mins(self):
        return self._mins

    @property
    def maxes(self):
        return self._maxes

    @property
    def lengths(self):
        return self._lengths

    @mins.setter
    def mins(self, mins):
        assert mins.shape == (3, )
        self._mins = mins
        self._lengths = self.maxes - self.mins

    @maxes.setter
    def maxes(self, maxes):
        assert maxes.shape == (3, )
        self._mins = maxes
        self._lengths = self.maxes - self.mins

    @lengths.setter
    def lengths(self, lengths):
        assert lengths.shape == (3, )
        self._maxes += lengths - self.lengths
        self._lengths = lengths

    def __repr__(self):
        return "Box(mins={}, maxes={})".format(self.mins, self.maxes)
