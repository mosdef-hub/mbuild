import numpy as np

from mbuild.coordinate_transform import translate, y_axis_transform
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.testing.tools import get_fn


class Initiator(Compound):
    def __init__(self):
        super(Initiator, self).__init__(self)

        # Look for data file in same directory as this python module.
        self.append_from_file(get_fn('initiator.pdb'))

        # Transform the coordinate system such that the two carbon atoms
        # that are part of the backbone are on the y axis, C_1 at the origin.
        y_axis_transform(self, new_origin=self.atoms[0], point_on_y_axis=self.atoms[21])

        # Add bottom port
        self.add(Port(anchor=self.atoms[0]), 'down')
        # Place the port.
        translate(self.down, self.atoms[0] + np.array([0.0, -0.07, 0.0]))

        # Add top port.
        self.add(Port(anchor=self.atoms[21]), 'up')
        # Place the port.
        translate(self.up, self.atoms[21] + np.array([0.0, 0.07, 0.0]))

if __name__ == "__main__":
    m = Initiator()


