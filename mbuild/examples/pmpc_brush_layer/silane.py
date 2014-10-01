import numpy as np

from mbuild.coordinate_transform import translate, x_axis_transform
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.testing.tools import get_fn


class Silane(Compound):

    def __init__(self, ):
        super(Silane, self).__init__()

        # Look for data file in same directory as this python module.
        self.append_from_file(get_fn('silane.pdb'))

        # Transform the coordinate system such that the silicon atom is at the origin
        # and the oxygen atoms are on the x axis.
        x_axis_transform(self, new_origin=self.SI[0], point_on_x_axis=self.O[0])

        # Add bottom port.
        self.add(Port(anchor=self.SI[0]), 'bottom_port')
        translate(self.bottom_port, np.array([0, -.07, 0]))

        # # Add top port.
        self.add(Port(anchor=self.SI[0]), 'top_port')
        translate(self.top_port, np.array([0, .07, 0]))

if __name__ == "__main__":
    m = Silane()
