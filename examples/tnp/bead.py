from mbuild.atom import Atom

__author__ = 'sallai'

from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import *

class Bead(Compound):
    """ """
    def __init__(self, particle_kind="particle"):
        Compound.__init__(self)

        self.add(Atom(pos=[0,0,0], kind=particle_kind), "particle")

        self.add(Port(anchor=self.particle), 'up')
        # rotate_around_z(self.up, np.pi)
        translate(self.up, np.array([0,0.7,0]))

        self.add(Port(anchor=self.particle), 'down')
        translate(self.down, np.array([0,-0.7,0]))

if __name__ == '__main__':
    m = Bead(particle_kind="bead")

    from mbuild.plot import Plot
    Plot(m, verbose=True, atoms=True, bonds=True, angles=False, dihedrals=False).show()


