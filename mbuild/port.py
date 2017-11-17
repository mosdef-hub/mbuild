import numpy as np

from mbuild.compound import Compound, Particle
from mbuild.coordinate_transform import unit_vector, angle
from mbuild import clone


class Port(Compound):
    """A set of four ghost Particles used to connect parts.

    Parameters
    ----------
    anchor : mb.Particle, optional, default=None
        A Particle associated with the port. Used to form bonds.
    orientation : array-like, shape=(3,), optional, default=[0, 1, 0]
        Vector along which to orient the port
    separation : float, optional, default=0
        Distance to shift port along the orientation vector from the anchor
        particle position. If no anchor is provided, the port will be shifted
        from the origin.

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
    def __init__(self, anchor=None, orientation=[0, 1, 0], separation=0):
        super(Port, self).__init__(name='Port', port_particle=True)
        self.anchor = anchor

        up = Compound(name='subport', port_particle=True)
        up.add(Particle(name='G', pos=[0.005, 0.0025, -0.0025],
                        port_particle=True), 'middle')
        up.add(Particle(name='G', pos=[0.005, 0.0225, -0.0025],
                        port_particle=True), 'top')
        up.add(Particle(name='G', pos=[-0.015, -0.0075, -0.0025],
                        port_particle=True), 'left')
        up.add(Particle(name='G', pos=[0.005, -0.0175, 0.0075],
                        port_particle=True), 'right')

        down = clone(up)
        down.rotate(np.pi, [0, 0, 1])

        '''
        Ports along antiparallel vectors should have 'up' and 'down' subports
        at equivalent positions, except with reverse labels (i.e. the 'up'
        subport for a port oriented in one direction will have the same xyz
        positions as the 'down' port oriented in the opposite direction. As
        a bit of a trick, when the z coordinate of the orientation vector is
        negative the subports are created as if the vector was of the opposite
        orientation and the labels are switched manually.
        '''
        swap_labels = False
        if orientation[2] < 0:
            swap_labels = True
            orientation = -np.array(orientation)
            self.add(up, 'down')
            self.add(down, 'up')
        else:
            self.add(up, 'up')
            self.add(down, 'down')
        self.used = False

        default_direction = [0, 1, 0]
        if np.allclose(
                np.asarray(default_direction), unit_vector(-np.asarray(orientation))):
            self.rotate(np.pi, [0, 0, 1])
        elif np.allclose(
                np.asarray(default_direction), unit_vector(np.asarray(orientation))):
            pass
        else:
            normal = np.cross(default_direction, orientation)
            self.rotate(angle(default_direction, orientation), normal)

        if anchor:
            self.translate_to(anchor.pos)

        if swap_labels == True:
            self.translate(separation*unit_vector(-orientation))
        else:
            self.translate(separation*unit_vector(orientation))

    def _clone(self, clone_of=None, root_container=None):
        newone = super(Port, self)._clone(clone_of, root_container)
        newone.anchor = clone(self.anchor, clone_of, root_container)
        newone.used = self.used
        return newone

    @property
    def center(self):
        """The cartesian center of the Port"""
        return np.mean(self.xyz_with_ports, axis=0)

    @property
    def direction(self):
        """The unit vector pointing in the 'direction' of the Port
        """
        return unit_vector(self.xyz_with_ports[1]-self.xyz_with_ports[0])
