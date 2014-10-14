import numpy as np

from mbuild.coordinate_transform import y_axis_transform, translate
from mbuild.compound import Compound
from mbuild.port import Port


class Cco(Compound):

    def __init__(self, ):
        super(Cco, self).__init__()

        self.append_from_file('cco.pdb')

        # Transform the coordinate system such that the silicon atom is at the origin
        # and the oxygen atoms are on the x axis.
        y_axis_transform(self, new_origin=self.C[0], point_on_y_axis=self.O[0])

        # Add bottom port.
        self.add(Port(anchor=self.O[0]), 'down')
        translate(self.down, self.O[0] + np.array([0, .07, 0]))

        # # Add top port.
        self.add(Port(anchor=self.C[0]), 'up')
        translate(self.up, self.C[0] - np.array([0, .07, 0]))

if __name__ == "__main__":
    m = Cco()
    m.visualize()
