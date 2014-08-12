__author__ = 'sallai'
import numpy as np

class Box(object):

    def __init__(self, lengths=None, mins=None, maxes=None):
        if lengths is not None:
            assert(mins is None and maxes is None)
            assert(np.shape(lengths) == (3,))
            self.mins = np.array([0,0,0])
            self.lengths = np.array(lengths)
            self.maxes = np.array(lengths)

        if maxes is not None:
            assert(mins is not None and lengths is None)
            assert(np.shape(mins) == (3,))
            assert(np.shape(maxes) == (3,))
            self.mins = np.array(mins)
            self.maxes = np.array(maxes)
            self.lengths = maxes - mins
