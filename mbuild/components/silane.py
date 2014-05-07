import sys
import os

import numpy as np
from numpy import pi

from mbuild.coordinate_transform import *
from mbuild.xyz import Xyz
from mbuild.compound import Compound
from mbuild.port import Port


class Silane(Compound):

    def __init__(self, ctx={}, alpha=0, ff='opls'):
        super(Silane, self).__init__(ctx=ctx)

        # Look for xyz file in same directory as this file.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        xyz_path = os.path.join(current_dir, 'silane.xyz')

        # Read xyz file with labels.
        # First C is labeled C_0, second C is C_1 etc.
        silane = Xyz(xyz_path, labels=True)
        self.add(silane, 'silane_xyz')

        # Pick two atoms to create axis on.
        if ff == 'opls':
            si_pos = silane.SI_0.pos
            o_pos = silane.opls_154_0.pos

        # Transform the coordinate system of silane such that the Si and
        # one oxygen are on the x axis, Si at the origin.
        silane.transform(AxisTransform(new_origin=si_pos,
                                       point_on_x_axis=o_pos))
        silane.transform(RotationAroundZ(pi / 2))
        if ff == 'opls':
            o_pos = silane.SI_0.pos
            si_pos = silane.opls_154_0.pos

        # Add bottom port.
        self.add(Port(), 'bottom_port')
        bottom_port_pos = si_pos - 1.5*(si_pos - o_pos)
        self.bottom_port.transform(Translation(bottom_port_pos))
        self.bottom_port.transform(RotationAroundX(pi / 2))

        # Add top port.
        self.add(Port(), 'top_port')
        top_port_pos = o_pos - 0.6*(o_pos - si_pos)
        self.top_port.transform(Translation(top_port_pos))
        self.top_port.transform(RotationAroundX(pi / 2))

if __name__ == "__main__":
    m = Silane()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()
