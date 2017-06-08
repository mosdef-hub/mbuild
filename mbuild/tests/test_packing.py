import pytest
import numpy as np

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.tests.base_test import BaseTest

class TestPacking(BaseTest):

    def test_fill_box(self, h2o):
        filled = mb.fill_box(h2o, n_compounds=50, box=[2, 2, 2, 4, 4, 4])
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2

    def test_fill_region(self, h2o):
        filled = mb.fill_region(h2o, n_compounds=50,
                                region=[3, 2, 2, 4, 4, 3])
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2
        assert np.min(filled.xyz[:,0]) >= 3
        assert np.max(filled.xyz[:,2]) <= 3

    def test_fill_region_multiple(self, ethane, h2o):
        filled = mb.fill_region(compound=[ethane, h2o], n_compounds=[2, 2],
                                region=[[2, 2, 2, 4, 4, 4], [4, 2, 2, 6, 4, 4]])
        assert filled.n_particles == 2 * 8 + 2 * 3
        assert filled.n_bonds == 2 * 7 + 2 * 2
        assert np.max(filled.xyz[:16, 0]) < 4
        assert np.min(filled.xyz[16:, 0]) > 4

    def test_fill_box_multiple(self, ethane, h2o):
        n_solvent = 100
        filled = mb.fill_box([ethane, h2o], [1, 100], box=[4, 4, 4])
        assert filled.n_particles == 8 + n_solvent * 3
        assert filled.n_bonds == 7 + n_solvent * 2
        assert len(filled.children) == 101

    def test_solvate(self, ethane, h2o):
        n_solvent = 100
        solvated = mb.solvate(ethane, h2o, n_solvent=n_solvent, box=[4, 4, 4])
        assert solvated.n_particles == 8 + n_solvent * 3
        assert solvated.n_bonds == 7 + n_solvent * 2

    def test_fill_box_seed(self, h2o):
        filled = mb.fill_box(h2o, n_compounds=50, box=[2, 2, 2])
        filled_same = mb.fill_box(h2o, n_compounds=50, box=[2, 2, 2])
        filled_diff = mb.fill_box(h2o, n_compounds=50, box=[2, 2, 2], seed=2)

    def test_solvate_multiple(self, methane, ethane, h2o):
        init_box = mb.fill_box(methane, 2, box=[4, 4, 4])
        solvated = mb.solvate(init_box, [ethane, h2o], [20, 20], box=[4, 4, 4])
        assert solvated.n_particles == 2*5 + 20*8 + 20*3
        assert len(solvated.children) == 41

    def test_fill_box_seed(self, ethane):
        filled = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2])
        filled_same = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2])
        filled_diff = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2], seed=2)
        assert np.array_equal(filled.xyz,filled_same.xyz)
        assert not np.array_equal(filled.xyz,filled_diff.xyz)

    def test_wrong_box(self, h2o):
        with pytest.raises(MBuildError):
            filled = mb.fill_box(h2o, n_compounds=50, box=[2, 2])
        with pytest.raises(MBuildError):
            filled = mb.fill_box(h2o, n_compounds=50, box=[2, 2, 2, 2])
