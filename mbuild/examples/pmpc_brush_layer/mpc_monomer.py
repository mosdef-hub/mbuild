from mbuild.coordinate_transform import translate, y_axis_transform, rotate_around_y
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.testing.tools import get_fn


class MpcMonomer(Compound):
    """A 2-(methacryloloxy) ethyl phosophorylcholine monomer."""
    def __init__(self, alpha=0):
        super(MpcMonomer, self).__init__(self)

        # Look for data file in same directory as this python module.
        self.append_from_file(get_fn('mpc.pdb'))

        # Transform the coordinate system of mpc such that the two carbon atoms
        # that are part of the backbone are on the y axis, c_backbone at the origin.
        C_top = self.atoms[37]
        # this can be achieved with the following alternative syntax:
        # C_top = self.labels["atom[37]"]
        # C_top = self.labels["atom"][37]
        C_bottom = self.atoms[1]

        y_axis_transform(self, new_origin=C_top, point_on_y_axis=C_bottom)

        # Add top port.
        self.add(Port(anchor=C_top), 'up')
        translate(self.up, C_top - (C_top - C_bottom)*1.50)

        # Add bottom port
        self.add(Port(anchor=C_bottom), 'down')
        rotate_around_y(self.down, alpha)
        translate(self.down, C_bottom - (C_bottom - C_top)*1.50)

if __name__ == "__main__":
    monomer = MpcMonomer().to_trajectory()

    monomer.save(filename='mpc.xyz')
