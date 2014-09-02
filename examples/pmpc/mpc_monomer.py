import os, sys
from mbuild.coordinate_transform import *
from mbuild.file_formats.mol2file import load_mol2
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.treeview import TreeView


class MpcMonomer(Compound):

    def __init__(self, alpha=0):
        Compound.__init__(self)

        # Look for data file in same directory as this python module.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        new_path = os.path.join(current_dir, 'mpc.mol2')
        load_mol2(new_path, part=self)

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
        # translate(top_port, c_backbone.pos - (c_backbone.pos - ch2_backbone.pos)*1.5)
        # translate(self.top_port, self.C_top - (self.C_top - self.C_bottom)*1.5)
        translate(self.top_port, C_top - (C_top - C_bottom)*1.5)

        # # Add bottom port
        self.add(Port(anchor=C_bottom), 'bottom_port')
        rotate_around_y(self.bottom_port, alpha)
        # translate(self.bottom_port, self.C_bottom - (self.C_bottom - self.C_top)*1.5)
        translate(self.bottom_port, C_bottom - (C_bottom - C_top)*1.5)

if __name__ == "__main__":
    m = MpcMonomer()

    # pdb.set_trace()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()

    # TreeView(m).show()
    from mbuild.file_formats.mol2file import write_mol2
    from mbuild.plugins.system import FlatCompound
    write_mol2(FlatCompound(compound=m),'mpc_monomer_out.mol2')