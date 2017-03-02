import numpy as np

import mbuild as mb


class Bead(mb.Compound):
    """A point particle with two ports pointing in opposite directions. """
    def __init__(self, particle_kind="bead"):
        """Initialize a Bead object.

        Args:
            particle_kind (str): Descriptive name for the Bead.
        """
        super(Bead, self).__init__()

        self.add(mb.Particle(name=particle_kind), particle_kind)

        self.add(mb.Port(anchor=self.labels[particle_kind]), 'up')
        self['up'].translate(np.array([0, 0.7, 0]))

        self.add(mb.Port(anchor=self.labels[particle_kind]), 'down')
        self['down'].translate(np.array([0, -0.7, 0]))

if __name__ == '__main__':
    bead = Bead(particle_kind="bead")
    print(bead)
