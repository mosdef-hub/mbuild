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

        Args:
            mins (np.ndarray of floats, optional): Minimum x, y, z
                coordinates.
            maxes (np.ndarray of floats, optional): Maximum x, y, z coordinates. 
            lengths (np.ndarray of floats, optional): Box length in x, y and z
                directions. 

        """
        # For clarities sake, we may want to throw a more descriptive exception
        # here re-iterating the points made in the doc string.
        if lengths is not None:
            assert(mins is None and maxes is None)
            assert(np.shape(lengths) == (3, ))
            self.mins = np.array([0.0, 0.0, 0.0])
            self.lengths = np.array(lengths)
            self.maxes = np.array(lengths)

        elif maxes is not None:
            assert(mins is not None and lengths is None)
            assert(np.shape(mins) == (3, ))
            assert(np.shape(maxes) == (3, ))
            self.mins = np.array(mins)
            self.maxes = np.array(maxes)
            self.lengths = self.maxes - self.mins

    def __repr__(self):
        return "Box(mins={}, maxes={})".format(self.mins, self.maxes)
