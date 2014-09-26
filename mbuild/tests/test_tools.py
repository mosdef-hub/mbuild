#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.tools` module. """
from mbuild.examples.bilayer.eceramidens import ECeramideNS
import numpy as np


class TestTools:

    def test_center_of_mass(self):
        e_ceramide_ns = ECeramideNS().to_trajectory()

        from mdtraj import compute_center_of_mass
        mdtraj_com = compute_center_of_mass(e_ceramide_ns)

        from mbuild.tools import compute_center_of_mass
        mbuild_com = compute_center_of_mass(e_ceramide_ns)

        assert mdtraj_com.all() == mbuild_com.all()

    def test_inertia_tensor(self):
        e_ceramide_ns = ECeramideNS().to_trajectory()

        from mbuild.tools import compute_inertia_tensor

        I = compute_inertia_tensor(e_ceramide_ns)
        ref = np.array([[[ 0.50559790, -0.01752252,  0.01581463],
                         [-0.01752252,  0.47113379,  0.00253122],
                         [ 0.01581463,  0.00253122,  0.05501593]]])
        assert I.all() == ref.all()

if __name__ == "__main__":
    pass

