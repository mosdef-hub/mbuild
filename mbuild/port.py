import itertools

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
    def __init__(self, anchor=None, orientation=None, separation=0):
        super(Port, self).__init__(name='Port', port_particle=True)
        self.anchor = anchor

        default_direction = np.array([0, 1, 0])
        if orientation is None:
            orientation = [0, 1, 0]
        orientation = np.asarray(orientation)

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
        self.add(up, 'up')
        self.add(down, 'down')
        self.used = False

        if np.allclose(
                default_direction, unit_vector(-orientation)):
            down.rotate(np.pi, [0, 0, 1])
            self.rotate(np.pi, [0, 0, 1])
        elif np.allclose(
                default_direction, unit_vector(orientation)):
            down.rotate(np.pi, [0, 0, 1])
        else:
            normal = np.cross(default_direction, orientation)
            self.rotate(angle(default_direction, orientation), normal)
            down.rotate(np.pi, normal)

        if anchor:
            self.translate_to(anchor.pos)

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

    @property
    def access_labels(self):
        """List of labels used to access the Port

        Returns
        -------
        list of str
            Strings that can be used to access this Port relative to self.root
        """
        access_labels = []
        for referrer in self.referrers:
            referrer_labels = [key for key, val in self.root.labels.items()
                               if val == referrer]
            port_labels = [key for key, val in referrer.labels.items()
                           if val == self]
            if referrer is self.root:
                for label in port_labels:
                    access_labels.append("['{}']".format(label))
            for label in itertools.product(referrer_labels, port_labels):
                access_labels.append("['{}']".format("']['".join(label)))

        return access_labels

    def __repr__(self):
        descr = list('<')
        descr.append(self.name + ', ')

        if self.anchor:
            descr.append("anchor: '{}', ".format(self.anchor.name))
        else:
            descr.append('anchor: None, ')

        descr.append('labels: {}, '.format(', '.join(self.access_labels)))

        descr.append('id: {}>'.format(id(self)))
        return ''.join(descr)
