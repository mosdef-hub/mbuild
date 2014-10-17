from numpy import sqrt, pi

from mbuild.compound import Compound
from mbuild.tools.polymer import Polymer
from mbuild.examples.tnp.sphere import Sphere
from mbuild.examples.tnp.bead import Bead
from mbuild.tools.mask import apply_mask, sphere_mask


class Tnp(Compound):
    """A spherical nanoparticle with tethered chains. """
    def __init__(self, ball_radius=10, n_chains=4, chain_length=10, monomer=None):
        """Initialize a tethered nanoparticle.

        Args:
            ball_radius (float): Radius of the nanoparticle.
            n_chains (int): Number of chains to attach to the nanoparticle.
            chain_length (int): Length of the chains being attached.
            monomer (Compound, optional): Type of chain being attached.
        """
        super(Tnp, self).__init__(self)

        if not monomer:
            monomer = Bead(particle_kind='t')

        n = 129  # TODO: make this tweakable
        self.add(Sphere(n=n, radius=ball_radius, port_distance_from_surface=0.7), "np")

        # generate 65 points on the surface of a unit sphere
        mask = sphere_mask(n_chains)

        # magnify it a bit
        mask *= ball_radius

        chain_proto = Polymer(monomer, n=chain_length)

        # apply chains to mask
        apply_mask(self.np, chain_proto , mask, guest_port_name="down")

        self.add_bond('np', 'np', sqrt(4 * ball_radius**2 * pi / n) - 0.5,
                                  sqrt(4 * ball_radius**2 * pi / n) + 0.5)
        self.add_bond('np', 't', 0.1, 0.3)
        self.add_bond('t', 'np', 0.1, 0.3)


def main():
    nano_particle = Tnp(n_chains=5, chain_length=10)
    nano_particle = nano_particle.to_trajectory()
    nano_particle.top.find_forcefield_terms()
    nano_particle.save(filename='tnp.hoomdxml')

if __name__ == "__main__":
    main()
