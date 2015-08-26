from copy import deepcopy

import numpy as np

from mbuild.atom import Atom
from mbuild.atom import Part
from mbuild.compound import Compound
from mbuild.coordinate_transform import rotate_around_z


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

        up = Compound(kind='subport')
        up.add(Atom(name='G', pos=[0, 0, 0]), 'middle')
        up.add(Atom(name='G', pos=[0, 0.02, 0]), 'top')
        up.add(Atom(name='G', pos=[-0.02, -0.01, 0]), 'left')
        up.add(Atom(name='G', pos=[0.0, -0.02, 0.01]), 'right')

        # down = deepcopy(up)
        down = up.clone()

        rotate_around_z(down, np.pi)

        self.add(up, 'up')
        self.add(down, 'down')

    def clone(self, clone_of=None, root_container=None, use_deepcopy=Part.USE_DEEPCOPY):
        if use_deepcopy:
            newone = super(Port, self).__deepcopy__(memo)
            newone.anchor = deepcopy(self.anchor, memo)
        else:
            if not clone_of:
                clone_of = dict()
            if not root_container:
                root_container = self
            newone = super(Port, self).clone(root_container=root_container, clone_of=clone_of, use_deepcopy=use_deepcopy)
            newone.anchor = self.anchor.clone(root_container=root_container, clone_of=clone_of, use_deepcopy=use_deepcopy)

        return newone

    def __deepcopy__(self, memo):
        newone = super(Port, self).__deepcopy__(memo)
        newone.anchor = deepcopy(self.anchor, memo)
        return newone
