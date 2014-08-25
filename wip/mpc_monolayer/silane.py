__author__ = 'sallai'

from mbuild.coordinate_transform import RotationAroundZ
from mbuild.file_formats.xyz import *
from mbuild.port import *


class Silane(Compound):

    def __init__(self, ctx={}, alpha=0):
        super(Silane, self).__init__(ctx=ctx)

        # read xyz file, turn on labeling (first C is labeled C_0, second C is C_1, and so on)
        silane = Xyz('silane.xyz', labels=True)
        #silane.transform(RotationAroundZ(pi/3))
        #silane.transform(Translation(np.array([2,3,4])))

        self.add(silane,'silane_xyz')

        # find the two atoms of the carbon chain
        si_pos = np.array(silane.Si_0.pos)
        o_pos = np.array(silane.O_0.pos)

        # transform the coordinate system of silane such that the Si and one oxygen are on the
        # the x axis, Si at the origin
        silane.transform(AxisTransform(new_origin=si_pos, point_on_x_axis=o_pos))
        silane.transform(RotationAroundZ(pi/2))

        o_pos = np.array(silane.Si_0.pos)
        si_pos = np.array(silane.O_0.pos)

        # add bottom port
        self.add(Port(), 'bottom_port')
        bottom_port_pos = si_pos - (si_pos - o_pos)*1.5
        self.bottom_port.transform(Translation(bottom_port_pos))
        self.bottom_port.transform(RotationAroundX(pi/2))

        # add top port
        self.add(Port(), 'top_port')
        top_port_pos = o_pos - (o_pos - si_pos)*0.7
        #self.top_port.transform(RotationAroundY(alpha))
        self.top_port.transform(Translation(top_port_pos))
        self.top_port.transform(RotationAroundX(pi/2))
if __name__ == "__main__":
    m = Silane()
    Plot(m, verbose=True).show()

    # TreeView(m).show()
