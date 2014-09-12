from math import atan2, asin

from mbuild.atom import Atom
from mbuild.compound import Compound
from mbuild.coordinate_transform import *
from mbuild.port import Port
from mbuild.plugins.mask import *

class Ball(Compound):

    def __init__(self, n=65, radius=1, port_distance_from_surface=.07):
        Compound.__init__(self)

        # generate 65 points on the surface of a unit sphere
        mask = sphere_mask(n)

        # magnify it a bit
        mask = mask * radius

        # create particles and ports at mask positions
        for i,pos in enumerate(mask):
            particle = Atom(kind="Si-cluster", pos=pos)
            self.add(particle,"Si-cluster_{}".format(i))
            port = Port(anchor=particle)
            self.add(port,"port_{}".format(i))

            # make the top of the port point toward the positive x axis
            rotate_around_z(port, -np.pi/2)
            # lift up (or down) the top of the port in the z direction
            rotate_around_y(port, -asin(pos[2]/radius))
            # rotate along the z axis
            rotate_around_z(port, atan2(pos[1], pos[0]))

            translate(port, pos + (pos/radius * port_distance_from_surface))

if __name__ == "__main__":
    m = Ball(n=65, radius=2)

    # pdb.set_trace()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()
