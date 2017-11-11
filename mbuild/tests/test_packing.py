import os

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

    def test_fill_box_density_box(self, h2o):
        filled = mb.fill_box(h2o, n_compounds=1000, density=1000)
        assert [3.1042931 < period < 3.1042932 for period in filled.periodicity]

    def test_fill_box_aspect_ratio(self, h2o):
        filled = mb.fill_box(h2o, n_compounds=1000,
                density=1000, aspect_ratio=[1, 2, 1])
        assert filled.periodicity[0]/filled.periodicity[1] == 0.5
        assert filled.periodicity[1]/filled.periodicity[2] == 2

    def test_fill_box_density_n_compounds(self, h2o):
        filled = mb.fill_box(h2o, density=1000,
                             box=mb.Box([3.1042931, 3.1042931, 3.1042931]))
        assert filled.n_particles == 3000

    def test_fill_box_compound_ratio(self, h2o, ethane):
        filled = mb.fill_box(compound=[h2o, ethane], density=800,
                compound_ratio=[2, 1], box=[2, 2, 2, 4, 4, 4])
        n_ethane = len([c for c in filled.children if c.name == 'Ethane'])
        n_water = len([c for c in filled.children if c.name == 'H2O'])
        assert n_water / n_ethane == 2

    def test_fill_region(self, h2o):
        filled = mb.fill_region(h2o, n_compounds=50,
                                region=[3, 2, 2, 4, 4, 3])
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2
        assert np.min(filled.xyz[:,0]) >= 3
        assert np.max(filled.xyz[:,2]) <= 3

    def test_fill_region_box(self, h2o):
        mybox = mb.Box([4, 4, 4])
        filled = mb.fill_region(h2o, n_compounds=50, region=mybox)
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2
        assert np.min(filled.xyz[:,0]) >= 0
        assert np.max(filled.xyz[:,2]) <= 4

    def test_fill_region_multiple(self, ethane, h2o):
        filled = mb.fill_region(compound=[ethane, h2o], n_compounds=[2, 2],
                                region=[[2, 2, 2, 4, 4, 4], [4, 2, 2, 6, 4, 4]])
        assert filled.n_particles == 2 * 8 + 2 * 3
        assert filled.n_bonds == 2 * 7 + 2 * 2
        assert np.max(filled.xyz[:16, 0]) < 4
        assert np.min(filled.xyz[16:, 0]) > 4

    def test_fill_region_multiple_boxes(self, ethane, h2o):
        box1 = mb.Box(mins=[2, 2, 2], maxs=[4, 4, 4])
        box2 = mb.Box(mins=[4, 2, 2], maxs=[6, 4, 4])
        filled = mb.fill_region(compound=[ethane, h2o], n_compounds=[2, 2],
                                region=[box1, box2])
        assert filled.n_particles == 2 * 8 + 2 * 3
        assert filled.n_bonds == 2 * 7 + 2 * 2
        assert np.max(filled.xyz[:16, 0]) < 4
        assert np.min(filled.xyz[16:, 0]) > 4

    def test_fill_region_multiple_types(self, ethane, h2o):
        box1 = mb.Box(mins=[2, 2, 2], maxs=[4, 4, 4])
        box2 = [4, 2, 2, 6, 4, 4]
        filled = mb.fill_region(compound=[ethane, h2o], n_compounds=[2, 2],
                                region=[box1, box2])
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

    def test_bad_args(self, h2o):
        with pytest.raises(ValueError):
            mb.fill_box(h2o, n_compounds=10)
        with pytest.raises(ValueError):
            mb.fill_box(h2o, density=1000)
        with pytest.raises(ValueError):
            mb.fill_box(h2o, box=[2, 2, 2])
        with pytest.raises(ValueError):
            mb.fill_box(h2o, n_compounds=10, density=1000, box=[2, 2, 2])
        with pytest.raises(ValueError):
            mb.fill_box(compound=[h2o, h2o], n_compounds=[10], density=1000)
        with pytest.raises(ValueError):
            mb.solvate(solute=h2o, solvent=[h2o], n_solvent=[10, 10], box=[2, 2, 2])
        with pytest.raises(ValueError):
            mb.fill_region(h2o, n_compounds=[10, 10], region=[2, 2, 2, 4, 4, 4])

    def test_write_temp_file(self, h2o):
        cwd = os.getcwd() # Must keep track of the temp dir that pytest creates
        filled = mb.fill_box(h2o, n_compounds=10, box=[4, 4, 4], temp_file='temp_file1.pdb')
        region = mb.fill_region(h2o, 10, [2, 2, 2, 4, 4, 4], temp_file='temp_file2.pdb')
        solvated = mb.solvate(filled, h2o, 10, box=[4, 4, 4], temp_file='temp_file3.pdb')
        assert os.path.isfile(os.path.join(cwd, 'temp_file1.pdb'))
        assert os.path.isfile(os.path.join(cwd, 'temp_file2.pdb'))
        assert os.path.isfile(os.path.join(cwd, 'temp_file3.pdb'))

    def test_packmol_error(self, h2o):
        with pytest.raises(RuntimeError):
            filled = mb.fill_box(h2o, n_compounds=10, box=[0, 0, 0])

    def test_packmol_warning(self, h2o):
        with pytest.warns(UserWarning):
            filled = mb.fill_box(h2o, n_compounds=10, box=[1, 1, 1], overlap=100)
