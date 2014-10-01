__author__ = 'sallai'

import numpy as np

from mbuild.atom import Atom
from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.coordinate_transform import translate


class Bead(Compound):
    """A point particle with two ports pointing in opposite directions. """
    def __init__(self, particle_kind="bead"):
        """Initialize a Bead object.

        Args:
            particle_kind (str): Descriptive name for the Bead.
        """
        Compound.__init__(self)

        self.add(Atom(kind=particle_kind), particle_kind)

        self.add(Port(anchor=self.labels[particle_kind]), 'up')
        translate(self.up, np.array([0, 0.7, 0]))

        self.add(Port(anchor=self.labels[particle_kind]), 'down')
        translate(self.down, np.array([0, -0.7, 0]))

if __name__ == '__main__':
    bead = Bead(particle_kind="bead")

    bead.visualize()

