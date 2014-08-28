import os, sys
from mbuild.file_formats.mol2file import load_mol2
from mbuild.coordinate_transform import *
from mbuild.compound import Compound
from mbuild.port import Port


class Silane(Compound):

    def __init__(self, ):
        super(Silane, self).__init__()

        # Look for data file in same directory as this python module.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        new_path = os.path.join(current_dir, 'silane.mol2')
        load_mol2(new_path, part=self)

        # Transform the coordinate system such that the silicon atom is at the origin
        # and the oxygen atoms are on the x axis.
        x_axis_transform(self, new_origin=self.SI[0], point_on_x_axis=self.O[0])

        # Add bottom port.
        self.add(Port(anchor=self.SI[0]), 'bottom_port')
        translate(self.bottom_port, np.array([0, -.7, 0]))

        # # Add top port.
        self.add(Port(anchor=self.SI[0]), 'top_port')
        translate(self.top_port, np.array([0, .7, 0]))

if __name__ == "__main__":
    m = Silane()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()
