from mbuild.coordinate_transform import *
from mbuild.file_formats.mol2file import load_mol2
from mbuild.compound import Compound
from mbuild.port import Port


class Initiator(Compound):

    def __init__(self):
        Compound.__init__(self)

        # Look for data file in same directory as this python module.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        new_path = os.path.join(current_dir, 'initiator.mol2')
        load_mol2(new_path, component=self)

        # Transform the coordinate system such that the two carbon atoms
        # that are part of the backbone are on the y axis, C_1 at the origin.
        y_axis_transform(self, new_origin=self.C_1, point_on_y_axis=self.C_22)

        # Add bottom port
        self.add(Port(anchor=self.C_1), 'bottom_port')
        # Place the port.
        translate(self.bottom_port, self.C_1 + np.array([0.0, -0.7, 0.0]))

        # Add top port.
        self.add(Port(anchor=self.C_22), 'top_port')
        # Place the port.
        translate(self.top_port, self.C_22 + np.array([0.0, 0.7, 0.0]))

if __name__ == "__main__":
    m = Initiator()
    from mbuild.plot import Plot
    # TreeView(m).show()
    Plot(m, verbose=True).show()
