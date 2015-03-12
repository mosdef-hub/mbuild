from copy import deepcopy

import numpy as np

from mbuild.atom import Atom
from mbuild.compound import Compound


class Port(Compound):
    """A set of four ghost Atoms used to connect parts.

    Parameters
    ----------
    anchor : mb.Atom, optional, default=None
        An atom associated with the port. Used to form bonds.

    Attributes
    ----------
    anchor : mb.Atom, optional, default=None
        An atom associated with the port. Used to form bonds.
    up : mb.Compound
        Collection of 4 ghost particles used to perform equivalence transforms.
        Faces the opposite direction as self.down.
    down : mb.Compound
        Collection of 4 ghost particles used to perform equivalence transforms.
        Faces the opposite direction as self.up.

    """
    def __init__(self, anchor=None):
        super(Port, self).__init__(kind='Port')
        self.anchor = anchor

        up = Compound()
        up.add(Atom(kind='G', pos=[0, 0, 0]), 'middle')
        up.add(Atom(kind='G', pos=[0, 0.02, 0]), 'top')
        up.add(Atom(kind='G', pos=[-0.02, -0.01, 0]), 'left')
        up.add(Atom(kind='G', pos=[0.0, -0.02, 0.01]), 'right')

        down = deepcopy(up)

        from mbuild.coordinate_transform import rotate_around_z
        rotate_around_z(down, np.pi)

        self.add(up, 'up')
        self.add(down, 'down')

    def __deepcopy__(self, memo):
        newone = super(Port, self).__deepcopy__(memo)
        newone.anchor = deepcopy(self.anchor, memo)
        return newone