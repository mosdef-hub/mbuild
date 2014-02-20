from mayavi.tools.sources import vector_scatter
from mbuild.coordinate_transform import RotationAroundY, RotationAroundZ
from mbuild.plot import Plot

__author__ = 'sallai'
from mbuild.compound import *
from mbuild.xyz import *
from mbuild.port import *

class MpcMonomer(Compound):

    def __init__(self, ctx={}):
        super(MpcMonomer, self).__init__(ctx=ctx)

        # read xyz file, turn on labeling (first C is labeled C_0, second C is C_1, and so on)
        mpc = Xyz('mpc_monomer.xyz', labels=True)
        self.add(mpc,'mpc_monomer_xyz')

        # find the two atoms of the carbon chain
        cbottom_pos = np.hstack(mpc.C_1.pos)
        ctop_pos = np.hstack(mpc.C_10.pos)

        print "ctop_pos=" +str(ctop_pos)
        print "cbottom_pos=" +str(cbottom_pos)

        # transform the coordinate system of mpc such that the two carbon atoms that are part of the backbone are on
        # the x axis, cbottom at the origin
        mpc.transform(AxisTransform(new_origin=cbottom_pos, point_on_x_axis=ctop_pos))

        # rotate MPC such that ctop and cbottom are on the y axis
        mpc.transform(RotationAroundZ(pi/2))


        # find the new positions of the two atoms of the carbon chain
        cbottom_pos = np.hstack(mpc.C_1.pos)
        ctop_pos = np.hstack(mpc.C_10.pos)

        print "ctop_pos=" +str(ctop_pos)
        print "cbottom_pos=" +str(cbottom_pos)


        # add bottom port
        self.add(Port(),'bottom_port')
        # compute position of bottom port
        bottom_port_pos = ctop_pos - (ctop_pos - cbottom_pos)*1.5
        # move the port there
        self.bottom_port.transform(Translation(bottom_port_pos))

        # add top port
        self.add(Port(),'top_port')
        # compute position of bottom port
        top_port_pos = cbottom_pos - (cbottom_pos - ctop_pos)*1.5
        # move the port there
        self.top_port.transform(Translation(top_port_pos))


if __name__ == "__main__":
    m = MpcMonomer()
    print m
    Plot(m, verbose=True).show()

    # TreeView(m).show()
