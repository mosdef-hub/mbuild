from __future__ import division

import numpy as np

from mbuild.trajectory import Trajectory
from mbuild.testing.tools import get_fn
from base_test import BaseTest

class TestCompute(BaseTest):

    def test_center_of_mass(self):
        e_ceramide_ns = Trajectory.load(get_fn('e-ceramide-ns.pdb'))

        import mdtraj as md
        mdtraj_com = md.compute_center_of_mass(e_ceramide_ns)
        from mbuild.tools.compute import compute_center_of_mass
        mbuild_com = compute_center_of_mass(e_ceramide_ns)
        assert (mdtraj_com == mbuild_com).all()

