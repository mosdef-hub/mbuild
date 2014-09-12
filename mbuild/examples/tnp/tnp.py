import numpy as np

from mbuild.compound import Compound
from mbuild.polymer import Polymer
from mbuild.examples.tnp.ball import Ball
from mbuild.examples.tnp.bead import Bead
from mbuild.plugins.mask import apply_mask, sphere_mask
from mbuild.tools import add_bond


class Tnp(Compound):

    def __init__(self, ball_radius=10, n_chains=4, chain_length=10, monomer=None):
        Compound.__init__(self)

        if not monomer:
            monomer = Bead()

        n = 129
        self.add(Ball(n=n, radius=ball_radius, port_distance_from_surface=0.7), "ball")

        # generate 65 points on the surface of a unit sphere
        mask = sphere_mask(n_chains)

        # magnify it a bit
        mask *= ball_radius

        chain_proto = Polymer(monomer, n=chain_length)

        # apply chains to mask
        apply_mask(self.ball, chain_proto , mask, guest_port_name="down")

        add_bond(self, "Si-cluster", "Si-cluster", np.sqrt(4*ball_radius*ball_radius*np.pi/n)-.5, np.sqrt(4*ball_radius*ball_radius*np.pi/n)+.5)
        add_bond(self, "Si-cluster", "particle", 0.1, 0.3)
        add_bond(self, "particle", "particle", 0.1, 0.3)


if __name__ == "__main__":
    from mbuild.formats.hoomdxml import save_hoomdxml
    nano_particle = Tnp(n_chains=5, chain_length=10)

    save_hoomdxml(nano_particle.to_trajectory(), filename='tnp.hoomdxml')

    #from mbuild.plot import Plot
    #Plot(nano_particle, verbose=False).show()
