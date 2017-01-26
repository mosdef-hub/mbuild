import pytest
import numpy as np

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.tests.base_test import BaseTest

class TestPacking(BaseTest):

    def test_fill_box(self, ethane):
        filled = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2, 4, 4, 4])
        assert filled.n_particles == 20 * 8
        assert filled.n_bonds == 20 * 7

    def test_solvate(self, ethane, h2o):
        n_solvent = 100
        solvated = mb.solvate(ethane, h2o, n_solvent=n_solvent, box=[4, 4, 4])
        assert solvated.n_particles == 8 + n_solvent * 3
        assert solvated.n_bonds == 7 + n_solvent * 2

    def test_fill_box_seed(self, ethane):
        filled = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2])
        filled_same = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2])
        filled_diff = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2], seed=2)
        assert np.array_equal(filled.xyz,filled_same.xyz)
        assert not np.array_equal(filled.xyz,filled_diff.xyz)

    def test_wrong_box(self, ethane):
        with pytest.raises(MBuildError):
            filled = mb.fill_box(ethane, n_compounds=20, box=[2, 2])
        with pytest.raises(MBuildError):
            filled = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2, 2])
