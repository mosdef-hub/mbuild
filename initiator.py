from mayavi.tools.sources import vector_scatter
from mbuild.coordinate_transform import RotationAroundY, RotationAroundZ
from mbuild.plot import Plot
from methane import Methane

__author__ = 'sallai'
from mbuild.compound import *
from mbuild.xyz import *
from mbuild.port import *
import pdb

class Initiator(Compound):

    def __init__(self, ctx={}, alpha=0):
        super(Initiator, self).__init__(ctx=ctx)

        # read xyz file, turn on labeling (first C is labeled C_0, second C is C_1, and so on)
        initiator = Xyz('initiator.xyz', labels=True)
        #initiator.transform(RotationAroundZ(pi/3))
        #initiator.transform(Translation(np.array([2,3,4])))

        self.add(initiator,'initiator_xyz')

        # find the two atoms of the carbon chain
        ctop_pos = np.array(initiator.C_7.pos)
        cbot_pos = np.array(initiator.C_0.pos)

        initiator.transform(AxisTransform(new_origin=cbot_pos, point_on_x_axis=ctop_pos))
        initiator.transform(RotationAroundZ(pi/2))

        ctop_pos = np.array(initiator.C_7.pos)
        cbot_pos = np.array(initiator.C_0.pos)

        # add bottom port
        self.add(Port(),'bottom_port')
        # compute position of bottom port
        bottom_port_pos = ctop_pos - (ctop_pos - cbot_pos)*1.1
        # move the port there
        self.bottom_port.transform(Translation(bottom_port_pos))

        # add top port
        self.add(Port(),'top_port')
        # rotate around y by alpha
        self.top_port.transform(RotationAroundY(alpha))

        # compute position of bottom port
        top_port_pos = cbot_pos - (cbot_pos - ctop_pos)*1.1
        # move the port there
        self.top_port.transform(Translation(top_port_pos))




if __name__ == "__main__":
    m = Initiator()
    print m
    Plot(m, verbose=True).show()

    # TreeView(m).show()
