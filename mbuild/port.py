from copy import deepcopy

import numpy as np

from mbuild.atom import Atom
from mbuild.compound import Compound
from mbuild.coordinate_transform import rotate_around_z


class Port(Compound):
    """A set of four ghost Atoms used to connect parts. """
    def __init__(self, anchor=None):
        super(Port, self).__init__(kind='Port')
        self.anchor = anchor

        up = Compound()
        up.add(Atom(kind='G', pos=np.array([0, 0, 0])), 'middle')
        up.add(Atom(kind='G', pos=np.array([0, 0.02, 0])), 'top')
        up.add(Atom(kind='G', pos=np.array([-0.02, -0.01, 0])), 'left')
        up.add(Atom(kind='G', pos=np.array([0.0, -0.02, 0.01])), 'right')

        down = deepcopy(up)
        rotate_around_z(down, np.pi)
        
        self.add(up, 'up')
        self.add(down, 'down')

    def __deepcopy__(self, memo):
        newone = super(Port, self).__deepcopy__(memo)
        newone.anchor = deepcopy(self.anchor, memo)
        return newone
