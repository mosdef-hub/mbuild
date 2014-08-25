from mbuild.coordinate_transform import RotationAroundY, RotationAroundZ

__author__ = 'sallai'
from mbuild.file_formats.xyz import *
from mbuild.port import *
import sys
import os


class MpcMonomer(Compound):

    def __init__(self, ctx={}, alpha=0):
        super(MpcMonomer, self).__init__(ctx=ctx)

        # read xyz file, turn on labeling (first C is labeled C_0, second C is C_1, and so on)
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))

        mpc = Xyz(os.path.join(current_dir, 'mpc_monomer.xyz'), labels=True)
        mpc.transform(RotationAroundZ(pi/3))
        mpc.transform(Translation(np.array([2,3,4])))

        self.add(mpc, 'mpc_monomer_xyz')

        # names of the backbone's carbons in the
        # c_1 = mpc.C_1
        # c_10 = mpc.C_10
        c_1 = mpc.opls_139_0
        c_10 = mpc.opls_136_0

        # find the two atoms of the carbon chain
        ctop_pos = np.hstack(c_1.pos) # C_1
        cbottom_pos = np.hstack(c_10.pos) # C_10

        # transform the coordinate system of mpc such that the two carbon atoms that are part of the backbone are on
        # the x axis, ctop at the origin
        mpc.transform(AxisTransform(new_origin=ctop_pos, point_on_x_axis=cbottom_pos))

        # rotate MPC such that cbottom and ctop are on the y axis
        mpc.transform(RotationAroundZ(pi/2))


        # find the new positions of the two atoms of the carbon chain
        ctop_pos = np.hstack(c_1.pos)
        cbottom_pos = np.hstack(c_10.pos)

        # print "cbottom_pos=" +str(cbottom_pos)
        # print "ctop_pos=" +str(ctop_pos)

        # add top port
        self.add(Port(),'top_port')
        # compute position of top port
        top_port_pos = cbottom_pos - (cbottom_pos - ctop_pos)*1.5
        # move the port there
        self.top_port.transform(Translation(top_port_pos))

        # add bottom port
        self.add(Port(),'bottom_port')
        # rotate around y by alpha
        self.bottom_port.transform(RotationAroundY(alpha))

        # compute position of top port
        bottom_port_pos = ctop_pos - (ctop_pos - cbottom_pos)*1.5
        # move the port there
        self.bottom_port.transform(Translation(bottom_port_pos))


if __name__ == "__main__":
    m = MpcMonomer()
    print m
    Plot(m, verbose=True).show()

    # TreeView(m).show()
