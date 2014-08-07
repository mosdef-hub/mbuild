from copy import deepcopy

__author__ = 'sallai'
import numpy as np

from compound import Compound
from atom import Atom

class Port(Compound):
    """ """
    def __init__(self, anchor=None):
        super(Port, self).__init__(kind='Port')
        self.anchor = anchor
        self.add(Atom(kind='G', pos=np.array([0, 0, 0])), 'middle')
        self.add(Atom(kind='G', pos=np.array([0, 0.2, 0])), 'top')
        self.add(Atom(kind='G', pos=np.array([-0.2, -0.1, 0])), 'left')
        self.add(Atom(kind='G', pos=np.array([0.0, -0.2, 0.1])), 'right')

    def __deepcopy__(self, memo):
        newone = super(Port, self).__deepcopy__(memo)
        newone.anchor = deepcopy(self.anchor, memo)

        return newone
