__author__ = 'sallai'
import numpy as np

from compound import Compound
from atom import Atom

class Port(Compound):
    """ """
    def __init__(self):
        super(Port, self).__init__(kind='Port')
        self.add(Atom(kind='G', pos=np.array([0, 0, 0])), 'middle')
        self.add(Atom(kind='G', pos=np.array([0, 0.2, 0])), 'top')
        self.add(Atom(kind='G', pos=np.array([-0.2, -0.1, 0])), 'left')
        self.add(Atom(kind='G', pos=np.array([0.0, -0.2, 0.1])), 'right')
