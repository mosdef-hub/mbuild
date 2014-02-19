from mayavi.tools.sources import vector_scatter
from mbuild.coordinate_transform import RotationAroundY
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

        # add bottom port
        self.add(Port(),'bottom_port')
        # rotate it by -90 degrees first
        self.bottom_port.transform(RotationAroundZ(-pi/2))
        # compute position of bottom port
        bottom_port_pos = ctop_pos - (ctop_pos - cbottom_pos)*1.5
        # transform the port's coordinate system such that bottom_port_pos is the new origin
        # and the x axis will point toward cbottom_pos
        self.bottom_port.transform(AxisTransform(new_origin=bottom_port_pos, point_on_x_axis=cbottom_pos))

        # add top port
        self.add(Port(),'top_port')
        # rotate it by -90 degrees first
        self.top_port.transform(RotationAroundZ(-pi/2))
        # compute position of top port
        top_port_pos = cbottom_pos + (ctop_pos - cbottom_pos)*1.5
        # compute a point that specifies where the x axis of the port's coordinate system should be transformed to
        top_port_x_axis = cbottom_pos + (ctop_pos - cbottom_pos)*2
        # transform the port's coordinate system such that top_port_pos is the new origin
        # and the x axis will point toward top_port_x_axis
        self.top_port.transform(AxisTransform(new_origin=top_port_pos, point_on_x_axis=top_port_x_axis))


if __name__ == "__main__":
    m = MpcMonomer()
    print m
    Plot(m, verbose=True).show()

    # TreeView(m).show()
