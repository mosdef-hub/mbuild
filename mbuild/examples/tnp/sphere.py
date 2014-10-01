import numpy as np

from mbuild.atom import Atom
from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.coordinate_transform import rotate_around_z, rotate_around_y, translate
from mbuild.tools.mask import sphere_mask


class Sphere(Compound):
    """A spherical arangement of particles with Ports. """
    def __init__(self, n=65, radius=1, port_distance_from_surface=.07):
        """Initialize a Sphere object.

        Args:
            n (int): Number of points used to construct the Sphere.
            radius (float): Radius of the Sphere.
            port_distance_from_surface (float): Distance of Ports from Sphere.
        """
        Compound.__init__(self)

        # Generate 65 points on the surface of a unit sphere.
        mask = sphere_mask(n)

        # Magnify the unit sphere by the provided radius.
        mask *= radius

        # Create particles and Ports at mask positions.
        for i, pos in enumerate(mask):
            particle = Atom(kind="np", pos=pos)
            self.add(particle, "np_{}".format(i))
            port = Port(anchor=particle)
            self.add(port, "port_{}".format(i))

            # Make the top of the port point toward the positive x axis.
            rotate_around_z(port, -np.pi/2)
            # Raise up (or down) the top of the port in the z direction.
            rotate_around_y(port, -np.arcsin(pos[2]/radius))
            # Rotate the Port along the z axis.
            rotate_around_z(port, np.arctan2(pos[1], pos[0]))
            # Move the Port a bit away from the surface of the Sphere.
            translate(port, pos + (pos/radius * port_distance_from_surface))

if __name__ == "__main__":
    m = Sphere(n=65, radius=2)
    m.visualize()
