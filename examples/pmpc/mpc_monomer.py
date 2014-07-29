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

        # Transform the coordinate system of mpc such that the two carbon atoms
        # that are part of the backbone are on the y axis, c_backbone at the origin.
        c_backbone = self.C_2
        ch2_backbone = self.C_38

        y_axis_transform(self, new_origin=c_backbone, point_on_y_axis=ch2_backbone)

        # Add top port.
        top_port = Port()
        self.add(top_port, 'top_port')
        # translate(top_port, c_backbone.pos - (c_backbone.pos - ch2_backbone.pos)*1.5)
        translate(top_port, c_backbone - (c_backbone - ch2_backbone)*1.5)

        # Add bottom port
        bottom_port = Port()
        self.add(bottom_port, 'bottom_port')
        rotate_around_y(bottom_port, alpha)
        translate(bottom_port, ch2_backbone - (ch2_backbone - c_backbone)*1.5)


if __name__ == "__main__":
    m = MpcMonomer()

    # pdb.set_trace()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()
