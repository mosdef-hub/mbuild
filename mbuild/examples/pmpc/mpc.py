import mbuild as mb


class MPC(mb.Compound):
    """A 2-(methacryloloxy) ethyl phosophorylcholine monomer."""
    def __init__(self, alpha=0):
        super(MPC, self).__init__()

        # Look for data file in same directory as this python module.
        mb.load('mpc.pdb', compound=self, relative_to_module=self.__module__)

        # Transform the coordinate system of mpc such that the two carbon atoms
        # that are part of the backbone are on the y axis, c_backbone at the origin.
        C_top = self[37]
        # this can be achieved with the following alternative syntax:
        # C_top = self.labels["atom[37]"]
        # C_top = self.labels["atom"][37]
        C_bottom = self[1]

        mb.y_axis_transform(self, new_origin=C_top, point_on_y_axis=C_bottom)

        # Add top port.
        self.add(mb.Port(anchor=C_top), label='up')
        mb.translate(self['up'], (C_bottom.pos - C_top.pos)*1.50)

        # Add bottom port
        self.add(mb.Port(anchor=C_bottom), 'down')
        mb.spin_y(self['down'], alpha)
        mb.translate(self['down'], (C_top.pos - C_bottom.pos)*1.50)

if __name__ == "__main__":
    monomer = MPC()
    print(monomer)
