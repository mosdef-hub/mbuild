from copy import deepcopy

import numpy as np

from mbuild.atom import Atom
from mbuild.compound import Compound
from mbuild.coordinate_transform import rotate_around_z
from mbuild import clone


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
        super(Port, self).__init__(name='Port')
        self.anchor = anchor

        up = Compound(name='subport')
        up.add(Atom(name='G', pos=[0, 0, 0]), 'middle')
        up.add(Atom(name='G', pos=[0, 0.02, 0]), 'top')
        up.add(Atom(name='G', pos=[-0.02, -0.01, 0]), 'left')
        up.add(Atom(name='G', pos=[0.0, -0.02, 0.01]), 'right')

        down = clone(up)

        rotate_around_z(down, np.pi)

        self.add(up, 'up')
        self.add(down, 'down')

    def _clone(self, clone_of=None, root_container=None):
        # if root_container is None:
        #     root_container = self
        # if clone_of is None:
        #     clone_of = dict()
        newone = super(Port, self)._clone(clone_of, root_container)
        newone.anchor = clone(self.anchor, clone_of, root_container)

        return newone

    @property
    def center(self):
        """The cartesian center of the Port"""
        return np.mean(self.xyz_with_ports, axis=0)
