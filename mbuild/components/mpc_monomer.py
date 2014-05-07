import sys
import os
import pdb

import numpy as np
from numpy import pi

from mbuild.coordinate_transform import *
from mbuild.xyz import Xyz
from mbuild.compound import Compound
from mbuild.port import Port


class MpcMonomer(Compound):

    def __init__(self, ctx={}, alpha=0, ff='opls'):
        super(MpcMonomer, self).__init__(ctx=ctx)

        # Look for xyz file in same directory as this file.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        xyz_path = os.path.join(current_dir, 'mpc_monomer.xyz')

        # Read xyz file with labels.
        # First C is labeled C_0, second C is C_1 etc.
        mpc = Xyz(xyz_path, labels=True)
        # mpc.transform(RotationAroundZ(pi / 3))
        # mpc.transform(Translation(np.array([2, 3, 4])))

        self.add(mpc, 'mpc_monomer_xyz')

        # Find two atoms of the carbon chain backbone.
        if ff == 'opls':
            c_backbone = mpc.opls_139_0
            ch2_backbone = mpc.opls_136_0

        ctop_pos = c_backbone.pos
        cbottom_pos = ch2_backbone.pos

        # Transform the coordinate system of mpc such that the two carbon atoms
        # that are part of the backbone are on the x axis, ctop at the origin.
        mpc.transform(AxisTransform(new_origin=ctop_pos,
                                    point_on_x_axis=cbottom_pos))

        # Rotate mpc such that cbottom and ctop are on the y axis.
        mpc.transform(RotationAroundZ(pi / 2))

        # Find the new positions of the two atoms of the carbon chain.
        ctop_pos = c_backbone.pos
        cbottom_pos = ch2_backbone.pos

        # Add top port.
        self.add(Port(), 'top_port')
        top_port_pos = cbottom_pos - (cbottom_pos - ctop_pos)*1.5
        self.top_port.transform(Translation(top_port_pos))

        # Add bottom port
        self.add(Port(), 'bottom_port')
        self.bottom_port.transform(RotationAroundY(alpha))
        bottom_port_pos = ctop_pos - (ctop_pos - cbottom_pos)*1.5
        self.bottom_port.transform(Translation(bottom_port_pos))


if __name__ == "__main__":
    m = MpcMonomer()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()
