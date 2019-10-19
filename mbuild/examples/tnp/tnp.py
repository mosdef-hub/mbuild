# -*- coding: utf-8 -*-


# -- ==tnp== --
from numpy import sqrt, pi

import mbuild as mb

from mbuild.examples.tnp.bead import Bead
from mbuild.examples.tnp.sphere import Sphere


class Tnp(mb.Compound):
    """A spherical nanoparticle with tethered chains. """
    def __init__(self, ball_radius=10, n_chains=4, chain_length=10, monomer=None):
        """Initialize a tethered nanoparticle.

        Args:
            ball_radius (float): Radius of the nanoparticle.
            n_chains (int): Number of chains to attach to the nanoparticle.
            chain_length (int): Length of the chains being attached.
            monomer (Compound, optional): Type of chain being attached.
        """
        super(Tnp, self).__init__()

        if not monomer:
            monomer = Bead(particle_kind='t')

        n = 129  # TODO: make this tweakable
        self.add(Sphere(n=n, radius=ball_radius, port_distance_from_surface=0.7), label="np")

        # Generate 65 points on the surface of a unit sphere.
        pattern = mb.SpherePattern(n_chains)

        # Magnify it a bit.
        pattern.scale(ball_radius)

        chain_proto = mb.Polymer(monomer, n=chain_length)

        # Apply chains to pattern.
        chain_protos, empty_backfill = pattern.apply_to_compound(chain_proto,
                guest_port_name="down", host=self['np'])
        self.add(chain_protos)

        self.generate_bonds('np', 'np', sqrt(4 * ball_radius ** 2 * pi / n) - 0.5,
                            sqrt(4 * ball_radius**2 * pi / n) + 0.5)
        self.generate_bonds('np', 't', 0.1, 0.3)
        self.generate_bonds('t', 'np', 0.1, 0.3)

# -- ==tnp== --