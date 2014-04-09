import sys
import os

import numpy as np
from numpy import pi

from mbuild.coordinate_transform import *
from mbuild.xyz import Xyz
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.plot import Plot


class Initiator(Compound):

    def __init__(self, ctx={}, alpha=0, ff='opls'):
        super(Initiator, self).__init__(ctx=ctx)

        # Look for xyz file in same directory as this file.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        xyz_path = os.path.join(current_dir, 'initiator.xyz')

        # Read xyz file with labels.
        # First C is labeled C_0, second C is C_1 etc.
        initiator = Xyz(xyz_path, labels=True)
        self.add(initiator, 'initiator_xyz')

        # Find two atoms of the carbon chain backbone.
        if ff == 'opls':
            cbottom_pos = np.hstack(initiator.opls_136_0.pos)
            ctop_pos = np.hstack(initiator.opls_139_0.pos)

        # Transform the coordinate system such that the two carbon atoms
        # that are part of the backbone are on the x axis, ctop at the origin.
        initiator.transform(AxisTransform(new_origin=cbottom_pos,
                                          point_on_x_axis=ctop_pos))
        # Rotate such that cbottom and ctop are on the y axis.
        initiator.transform(RotationAroundZ(pi / 2))

        # Add bottom port
        self.add(Port(), 'bottom_port')
        # Rotate around y by alpha to give the molecule a little twist.
        self.bottom_port.transform(RotationAroundY(alpha))
        # Place the port.
        cbottom_pos = np.hstack(initiator.opls_136_0.pos)
        bottom_port_pos = cbottom_pos + (0.0, -0.7, 0.0)
        self.bottom_port.transform(Translation(bottom_port_pos))

        # Add top port.
        self.add(Port(), 'top_port')
        ctop_pos = np.hstack(initiator.opls_139_0.pos)
        top_port_pos = ctop_pos + (0.0, 0.7, 0.0)
        self.top_port.transform(Translation(top_port_pos))

if __name__ == "__main__":
    m = Initiator()
    Plot(m, verbose=True).show()
