"""Ports used to facilitate bond formation."""
import itertools
from warnings import warn

import numpy as np

from mbuild import clone
from mbuild.compound import Compound, Particle
from mbuild.coordinate_transform import angle, unit_vector


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
        super(Port, self).__init__(name="Port", port_particle=True)
        self.anchor = anchor

        default_direction = np.array([0, 1, 0])
        if orientation is None:
            orientation = [0, 1, 0]
        orientation = np.asarray(orientation)

        up = Compound(name="subport", port_particle=True)
        pos = [0.005, 0.0025, -0.0025]
        up.add(Particle(name="G", pos=pos, port_particle=True), "middle")
        pos = [0.005, 0.0225, -0.0025]
        up.add(Particle(name="G", pos=pos, port_particle=True), "top")
        pos = [-0.015, -0.0075, -0.0025]
        up.add(Particle(name="G", pos=pos, port_particle=True), "left")
        pos = [0.005, -0.0175, 0.0075]
        up.add(Particle(name="G", pos=pos, port_particle=True), "right")

        down = clone(up)
        self.add(up, "up")
        self.add(down, "down")
        self.used = False

        if np.allclose(default_direction, unit_vector(-orientation)):
            down.rotate(np.pi, [0, 0, 1])
            self.rotate(np.pi, [0, 0, 1])
        elif np.allclose(default_direction, unit_vector(orientation)):
            down.rotate(np.pi, [0, 0, 1])
        else:
            normal = np.cross(default_direction, orientation)
            self.rotate(angle(default_direction, orientation), normal)
            down.rotate(np.pi, normal)

        if anchor:
            self.translate_to(anchor.pos)

        self.translate(separation * unit_vector(orientation))

    def _clone(self, clone_of=None, root_container=None):
        newone = super(Port, self)._clone(clone_of, root_container)
        newone.anchor = clone(self.anchor, clone_of, root_container)
        newone.used = self.used
        return newone

    def update_separation(self, separation):
        """Change the distance between a port and its anchor particle.

        separation : float, required
            Distance to shift port along the orientation vector from the anchor
            particle position. If no anchor is provided, the port will be
            shifted from the origin.
        """
        if self.used:
            warn(
                "This port is already being used and changing its separation "
                "will have no effect on the distance between particles."
            )

        if self.anchor:
            self.translate_to(self.anchor.pos)
        else:
            self.translate_to((0, 0, 0))

        self.translate(separation * self.direction)

    def update_orientation(self, orientation):
        """Change the direction between a port and its anchor particle.

        orientation : array-like, shape=(3,), required
            Vector along which to orient the port
        """
        if self.used:
            warn(
                "This port is already being used and changing its orientation "
                "will have no effect on the direction between particles."
            )

        orientation = np.asarray(orientation).reshape(3)
        init_separation = self.separation
        normal = np.cross(self.direction, orientation)
        # Move to origin to perform rotation
        self.translate_to((0, 0, 0))
        self.rotate(angle(self.direction, orientation), normal)
        self.labels["down"].rotate(np.pi, normal)
        self.labels["up"].rotate(np.pi, normal)
        # Move back to it's anchor particle
        self.update_separation(init_separation)

    @property
    def center(self):
        """Get the cartesian center of the Port."""
        return np.mean(self.xyz_with_ports, axis=0)

    @property
    def direction(self):
        """Get the unit vector pointing in the 'direction' of the Port."""
        return unit_vector(self.xyz_with_ports[1] - self.xyz_with_ports[0])

    @property
    def separation(self):
        """Get the distance between a port and its anchor particle.

        If the port has no anchor particle, returns None.
        """
        if self.anchor:
            return np.linalg.norm(self.center - self.anchor.pos)
        else:
            warn(
                "This port is not anchored to another particle. Returning a "
                "separation of None"
            )
            return None

    @property
    def access_labels(self):
        """List labels used to access the Port.

        Returns
        -------
        list of str
            Strings that can be used to access this Port relative to self.root
        """
        access_labels = []
        for referrer in self.referrers:
            referrer_labels = [
                key for key, val in self.root.labels.items() if val == referrer
            ]
            port_labels = [
                key for key, val in referrer.labels.items() if val == self
            ]
            if referrer is self.root:
                for label in port_labels:
                    access_labels.append(label)
            for label in itertools.product(referrer_labels, port_labels):
                access_labels.append(label)

        return access_labels

    def __repr__(self):
        """Return the Port's representation."""
        descr = list("<")
        descr.append(self.name + ", ")

        if self.anchor:
            descr.append("anchor: '{}', ".format(self.anchor.name))
        else:
            descr.append("anchor: None, ")

        descr.append("labels: {}, ".format(", ".join(self.access_labels)))

        descr.append("id: {}>".format(id(self)))
        return "".join(descr)
