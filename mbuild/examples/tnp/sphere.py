from numpy import pi, arctan2, arcsin

import mbuild as mb


class Sphere(mb.Compound):
    """A spherical arangement of particles with Ports. """
    def __init__(self, n=65, radius=1, port_distance_from_surface=.07):
        """Initialize a Sphere object.

        Args:
            n (int): Number of points used to construct the Sphere.
            radius (float): Radius of the Sphere.
            port_distance_from_surface (float): Distance of Ports from Sphere.
        """
        super(Sphere, self).__init__()
        particle = mb.Particle(name='np')
        particle.add(mb.Port(anchor=particle), label='out')

        # Generate 65 points on the surface of a unit sphere.
        pattern = mb.SpherePattern(n)
        # Magnify the unit sphere by the provided radius.
        pattern.scale(radius)

        particles = pattern.apply(particle, orientation='normal', compound_port='out')
        self.add(particles, label='np_[$]')

        # Create particles and Ports at pattern positions.
        for i, pos in enumerate(pattern.points):
            particle = mb.Particle(name="np", pos=pos)
            self.add(particle, "np_{}".format(i))
            port = mb.Port(anchor=particle)
            self.add(port, "port_{}".format(i))

            # Make the top of the port point toward the positive x axis.
            port.spin(-pi/2, [0, 0, 1])
            # Raise up (or down) the top of the port in the z direction.
            port.spin(-arcsin(pos[2]/radius), [0, 1, 0])
            # Rotate the Port along the z axis.
            port.spin(arctan2(pos[1], pos[0]), [0, 0, 1])
            # Move the Port a bit away from the surface of the Sphere.
            port.translate(pos/radius * port_distance_from_surface)

if __name__ == "__main__":
    sphere = Sphere(n=4, radius=2)
    print(sphere)
