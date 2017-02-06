import numpy as np

from mbuild.compound import Compound, Particle
from mbuild.coordinate_transform import rotate_around_z, translate_to
from mbuild import clone


class Port(Compound):
    """A set of four ghost Particles used to connect parts.

    Parameters
    ----------
    anchor : mb.Particle, optional, default=None
        A Particle associated with the port. Used to form bonds.

    Attributes
    ----------
    anchor : mb.Particle, optional, default=None
        A Particle associated with the port. Used to form bonds.
    up : mb.Compound
        Collection of 4 ghost particles used to perform equivalence transforms.
        Faces the opposite direction as self['down'].
    down : mb.Compound
        Collection of 4 ghost particles used to perform equivalence transforms.
        Faces the opposite direction as self['up'].
    used : bool
        Status of whether a port has been occupied following an equivalence
        transform.

    """
    def __init__(self, anchor=None):
        super(Port, self).__init__(name='Port', port_particle=True)
        self.anchor = anchor

        up = Compound(name='subport', port_particle=True)
        up.add(Particle(name='G', pos=[0, 0, 0], port_particle=True), 'middle')
        up.add(Particle(name='G', pos=[0, 0.02, 0], port_particle=True), 'top')
        up.add(Particle(name='G', pos=[-0.02, -0.01, 0], port_particle=True), 'left')
        up.add(Particle(name='G', pos=[0.0, -0.02, 0.01], port_particle=True), 'right')

        down = clone(up)
        rotate_around_z(down, np.pi)

        self.add(up, 'up')
        self.add(down, 'down')
        self.used = False

        if anchor:
            translate_to(self, anchor.pos)

    def _clone(self, clone_of=None, root_container=None):
        newone = super(Port, self)._clone(clone_of, root_container)
        newone.anchor = clone(self.anchor, clone_of, root_container)
        newone.used = self.used
        return newone

    @property
    def center(self):
        """The cartesian center of the Port"""
        return np.mean(self.xyz_with_ports, axis=0)
