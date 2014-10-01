from mbuild.coordinate_transform import translate, y_axis_transform, rotate_around_y
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.testing.tools import get_fn


class MpcMonomer(Compound):

    def __init__(self, alpha=0):
        Compound.__init__(self)

        # Look for data file in same directory as this python module.
        self.append_from_file(get_fn('mpc.pdb'))

        # Transform the coordinate system of mpc such that the two carbon atoms
        # that are part of the backbone are on the y axis, c_backbone at the origin.
        C_top = self.atom[37]
        # this can be achieved with the following alternative syntax:
        # C_top = self.labels["atom[37]"]
        # C_top = self.labels["atom"][37]
        C_bottom = self.atom[1]

        y_axis_transform(self, new_origin=C_top, point_on_y_axis=C_bottom)

        # Add top port.
        self.add(Port(anchor=C_top), 'top_port')
        translate(self.top_port, C_top - (C_top - C_bottom)*1.50)

        # Add bottom port
        self.add(Port(anchor=C_bottom), 'bottom_port')
        rotate_around_y(self.bottom_port, alpha)
        translate(self.bottom_port, C_bottom - (C_bottom - C_top)*1.50)

if __name__ == "__main__":
    monomer = MpcMonomer().to_trajectory()

    monomer.save(filename='mpc.xyz')
