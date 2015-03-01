#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.tools.compute` module. """
import pytest
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
        assert mdtraj_com.all() == mbuild_com.all()

    def test_inertia_tensor(self):
        e_ceramide_ns = Trajectory.load(get_fn('e-ceramide-ns.pdb'))

        from mbuild.tools.compute import compute_inertia_tensor
        I = compute_inertia_tensor(e_ceramide_ns)
        ref = np.array([[[ 0.50559790, -0.01752252,  0.01581463],
                         [-0.01752252,  0.47113379,  0.00253122],
                         [ 0.01581463,  0.00253122,  0.05501593]]])
        assert I.all() == ref.all()


