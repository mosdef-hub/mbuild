__author__ = 'sallai'

from copy import deepcopy

import numpy as np

from atom import Atom
from compound import Compound


class Port(Compound):
    """A set of four ghost Atoms used to connect parts. """
    def __init__(self, anchor=None):
        super(Port, self).__init__(kind='Port')
        self.anchor = anchor
        self.add(Atom(kind='G', pos=np.array([0, 0, 0])), 'middle')
        self.add(Atom(kind='G', pos=np.array([0, 0.02, 0])), 'top')
        self.add(Atom(kind='G', pos=np.array([-0.02, -0.01, 0])), 'left')
        self.add(Atom(kind='G', pos=np.array([0.0, -0.02, 0.01])), 'right')

    def __deepcopy__(self, memo):
        newone = super(Port, self).__deepcopy__(memo)
        newone.anchor = deepcopy(self.anchor, memo)
        return newone
