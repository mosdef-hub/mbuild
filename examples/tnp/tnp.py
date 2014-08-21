from math import atan2, asin
import sys
import os
import pdb

from numpy import pi
from examples.tnp.ball import Ball
from examples.tnp.bead import Bead

from mbuild.coordinate_transform import *
from mbuild.polymer import Polymer
from mbuild.xyz import load_xyz
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.tools import *

class Tnp(Compound):

    def __init__(self, ball_radius=10, n_chains=4, chain_length=10):
        Compound.__init__(self)

        self.add(Ball(n=129, radius=ball_radius, port_distance_from_surface=.7), "ball")

        # generate 65 points on the surface of a unit sphere
        mask = sphere_mask(n_chains)

        # magnify it a bit
        mask = mask * ball_radius

        chain_proto = Polymer(Bead(), n=chain_length)

        # apply chains to mask
        apply_mask(self.ball, chain_proto , mask, guest_port_name="down")


if __name__ == "__main__":
    m = Tnp(n_chains=5, chain_length=10)

    # pdb.set_trace()
    from mbuild.plot import Plot
    Plot(m, verbose=False).show()
