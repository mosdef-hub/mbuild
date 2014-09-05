from mbuild.examples.tnp.ball import Ball
from mbuild.examples.tnp.bead import Bead

from mbuild.polymer import Polymer
from mbuild.plugins.mask import *
from mbuild.tools import *

class Tnp(Compound):

    def __init__(self, ball_radius=10, n_chains=4, chain_length=10):
        Compound.__init__(self)

        n = 129
        self.add(Ball(n=n, radius=ball_radius, port_distance_from_surface=.7), "ball")

        # generate 65 points on the surface of a unit sphere
        mask = sphere_mask(n_chains)

        # magnify it a bit
        mask = mask * ball_radius

        chain_proto = Polymer(Bead(), n=chain_length)

        # apply chains to mask
        apply_mask(self.ball, chain_proto , mask, guest_port_name="down")

        add_bond(self, "particle", "particle", np.sqrt(4*ball_radius*ball_radius*np.pi/n)-.5, np.sqrt(4*ball_radius*ball_radius*np.pi/n)+.5)


if __name__ == "__main__":
    m = Tnp(n_chains=5, chain_length=10)

    # pdb.set_trace()
    from mbuild.plot import Plot
    Plot(m, verbose=False).show()
