import sys
import os
import pdb

from numpy import pi

from mbuild.coordinate_transform import *
from mbuild.mol2file import load_mol2
from mbuild.compound import Compound
from mbuild.port import Port

class MpcMonomer(Compound):

    def __init__(self, alpha=0):
        Compound.__init__(self)

        # Look for data file in same directory as this python module.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        new_path = os.path.join(current_dir, 'mpc.mol2')
        load_mol2(new_path, component=self)

        # print(self.labels)

        c_backbone = self.labels['C.3_2']
        ch2_backbone = self.labels['C.3_38']
        ctop_pos = c_backbone.pos
        cbottom_pos = ch2_backbone.pos

        # Transform the coordinate system of mpc such that the two carbon atoms
        # that are part of the backbone are on the x axis, ctop at the origin.
        transform(self, AxisTransform(new_origin=ctop_pos,
                                    point_on_x_axis=cbottom_pos))

        # Rotate mpc such that cbottom and ctop are on the y axis.
        transform(self, RotationAroundZ(pi / 2))

        # Find the new positions of the two atoms of the carbon chain.
        ctop_pos = c_backbone.pos
        cbottom_pos = ch2_backbone.pos

        # Add top port.
        top_port = Port()
        self.add(top_port, 'top_port')
        top_port_pos = cbottom_pos - (cbottom_pos - ctop_pos)*1.5
        transform(top_port, Translation(top_port_pos))

        # Add bottom port
        bottom_port = Port()
        self.add(bottom_port, 'bottom_port')
        transform(bottom_port, RotationAroundY(alpha))
        bottom_port_pos = ctop_pos - (ctop_pos - cbottom_pos)*1.5
        transform(bottom_port, Translation(bottom_port_pos))


if __name__ == "__main__":
    m = MpcMonomer()


    # pdb.set_trace()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()
