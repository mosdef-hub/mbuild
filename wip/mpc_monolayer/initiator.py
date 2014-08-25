from mbuild.coordinate_transform import RotationAroundY, RotationAroundZ

__author__ = 'sallai'
from mbuild.file_formats.xyz import *
from mbuild.port import *

class Initiator(Compound):

    def __init__(self, ctx={}, alpha=0):
        super(Initiator, self).__init__(ctx=ctx)

        # read xyz file, turn on labeling (first C is labeled C_0, second C is C_1, and so on)
        c = Xyz('initiator.xyz', labels=True)
        self.add(c,'initiator_xyz')

        # find two atoms of the carbon chain backbone
        cbottom_pos = np.hstack(c.C_7.pos)
        ctop_pos = np.hstack(c.C_6.pos)

        # transform the coordinate system of c such that the two carbon atoms that are part of the backbone are on
        # the x axis, ctop at the origin
        c.transform(AxisTransform(new_origin=cbottom_pos, point_on_x_axis=ctop_pos))
        # rotate MPC such that cbottom and ctop are on the y axis
        c.transform(RotationAroundZ(pi/2))

        # add bottom port
        self.add(Port(),'bottom_port')
        # rotate around y by alpha to give the molecule a little twist
        self.bottom_port.transform(RotationAroundY(alpha))
        # compute position of top port
        bottom_port_pos = (0.0, -0.7, 0.0)
        # # move the port there
        self.bottom_port.transform(Translation(bottom_port_pos))


        # find the position of the c atom at the top end of the backbone
        ctop_pos = np.hstack(c.C_0.pos)
        # add top port
        self.add(Port(),'top_port')
        # compute position of top port
        top_port_pos = ctop_pos + (0.0, 0.7, 0.0)
        # move the port there
        self.top_port.transform(Translation(top_port_pos))


if __name__ == "__main__":
    m = Initiator()
    print m
    Plot(m, verbose=True).show()
    # TreeView(m).show()
