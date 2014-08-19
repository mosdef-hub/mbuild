from math import atan2, asin
import sys
import os
import pdb

from numpy import pi
from examples.tnp.bead import Bead

from mbuild.coordinate_transform import *
from mbuild.polymer import Polymer
from mbuild.xyz import load_xyz
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.tools import *

class Ball(Compound):

    def __init__(self, n=65, radius=10):
        Compound.__init__(self)

        # generate 65 points on the surface of a unit sphere
        mask = sphere_mask(n)

        # magnify it a bit
        mask = mask * radius

        # create particles and ports at mask positions
        for i,pos in enumerate(mask):
            particle = Atom(kind="particle", pos=pos)
            self.add(particle,"particle_{}".format(i))
            port = Port(anchor=particle)
            self.add(port,"port_{}".format(i))
            rotate_around_z(port, np.pi)
            rotate_around_x(port, -asin(pos[2]/radius))
            rotate_around_z(port, +np.pi/2 + atan2(pos[1], pos[0]))
            translate(port, pos * 1.1)

if __name__ == "__main__":
    m = Ball(n=128, radius = 20)

    # pdb.set_trace()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()
