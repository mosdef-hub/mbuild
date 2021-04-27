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

    def test_fill_sphere(self, h2o):
        filled = mb.fill_sphere(h2o, sphere=[3, 3, 3, 1.5], n_compounds=50)
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2

        center = np.array([3.0, 3.0, 3.0])
        assert np.alltrue(np.linalg.norm(filled.xyz - center, axis=1) < 1.5)

    def test_fill_sphere_density(self, h2o):
        filled = mb.fill_sphere(h2o, sphere=[3, 3, 3, 1.5], density=1000)
        assert filled.n_particles == 921

    def test_fill_sphere_compound_ratio(self, h2o, ethane):
        filled = mb.fill_sphere(compound=[h2o, ethane], sphere=[3, 3, 3, 1.5],
                density=800, compound_ratio=[2, 1])
        n_ethane = len([c for c in filled.children if c.name == 'Ethane'])
        n_water = len([c for c in filled.children if c.name == 'H2O'])
        assert n_water / n_ethane == 2

    def test_fill_sphere_bad_args(self, h2o, ethane):
        with pytest.raises(ValueError):
            mb.fill_sphere(compound=h2o, sphere=[4, 4, 4, 1])
        with pytest.raises(ValueError):
            mb.fill_sphere(compound=h2o, n_compounds=100,
                density=100, sphere=[4, 4, 4, 1])
        with pytest.raises(ValueError):
            mb.fill_sphere(compound=h2o, density=1000, sphere='yes')
        with pytest.raises(ValueError):
            mb.fill_sphere(compound=[h2o, ethane], n_compounds=1000, sphere=[1, 1, 1, 4])
        with pytest.raises(ValueError):
            mb.fill_sphere(compound=h2o, n_compounds=[10, 10], sphere=[1, 1, 1, 4])
        with pytest.raises(ValueError):
            mb.fill_sphere(compound=h2o, n_compounds=100, sphere=[1, 1, 1, 4])


    def test_fill_region(self, h2o):
        filled = mb.fill_region(h2o, n_compounds=50,
                                region=[3, 2, 2, 5, 5, 5])
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2
        assert np.min(filled.xyz[:,0]) >= 3
        assert np.min(filled.xyz[:,1]) >= 2
        assert np.min(filled.xyz[:,2]) >= 2
        assert np.max(filled.xyz[:,0]) <= 5
        assert np.max(filled.xyz[:,1]) <= 5
        assert np.max(filled.xyz[:,2]) <= 5

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
        with pytest.raises(ValueError):
            mb.fill_box(compound=[h2o, h2o], n_compounds=[10], density=1000,
                        fix_orientation=[True, True, True])

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

    def test_packmol_log_error(self, h2o):
        try:
            filled = mb.fill_box(h2o, n_compounds=10, box=[0, 0, 0])
        except(RuntimeError):
            with open("log.txt", "r") as logfile:
                assert "ERROR" in logfile.read()

    def test_packmol_warning(self, h2o):
        with pytest.warns(UserWarning):
            filled = mb.fill_box(h2o, n_compounds=10, box=[1, 1, 1], overlap=10)

    def test_rotate(self, h2o):
        filled = mb.fill_box(h2o, 2, box=[1, 1, 1], fix_orientation=True)
        w0 = filled.xyz[:3]
        w1 = filled.xyz[3:]
        # Translate w0 and w1 to COM
        w0 -= w0.sum(0) / len(w0)
        w1 -= w1.sum(0) / len(w1)
        assert np.isclose(w0, w1).all()

    def test_no_rotate(self, h2o):
        filled = mb.fill_box([h2o, h2o], [1, 1], box=[1, 1, 1], fix_orientation=[False, True])
        w0 = filled.xyz[:3]
        w1 = filled.xyz[3:]
        # Translate w0 and w1 to COM
        w0 -= w0.sum(0) / len(w0)
        w1 -= w1.sum(0) / len(w1)
        assert np.isclose(w0, w1).all() is not True

    def test_remove_port(self):
        from mbuild.lib.recipes import Alkane

        butane = Alkane(n=4)
        butane.remove(butane[-1])
        box = mb.fill_box(butane, n_compounds=10, density=1)

    def test_sidemax(self):
        from mbuild.lib.molecules import Methane
        ch4 = Methane()
        #With default sidemax
        box_of_methane = mb.fill_box(ch4,
                    box=[1000, 1000, 1000],
                    n_compounds=500)
        sphere_of_methane = mb.fill_sphere(ch4,
                    sphere=[1000, 1000, 1000, 1000],
                    n_compounds=500)
        assert all(box_of_methane.boundingbox.lengths < [110, 110, 110])
        assert all(sphere_of_methane.boundingbox.lengths < [210, 210, 210])

        #With adjusted sidemax
        big_box_of_methane = mb.fill_box(ch4,
                    box=[1000, 1000, 1000],
                    n_compounds=500,
                    sidemax=1000.0)
        big_sphere_of_methane = mb.fill_sphere(ch4,
                    sphere=[1000, 1000, 1000, 1000],
                    n_compounds=500,
                    sidemax=2000.0)
        assert all(big_box_of_methane.boundingbox.lengths > [900, 900, 900])
        assert all(big_sphere_of_methane.boundingbox.lengths > [1800, 1800, 1800])

    def test_box_edge(self, h2o, methane):
        system_box = mb.Box(lengths=(1.8, 1.8, 1.8))
        packed = mb.fill_box(compound=h2o,
                             n_compounds=100,
                             box=system_box,
                             edge=0.2)
        edge_sizes = system_box.lengths - packed.boundingbox.lengths
        assert np.allclose(edge_sizes, np.array([0.4]*3), atol=0.1)

        region = mb.fill_region(compound=h2o,
                                n_compounds=100,
                                region=system_box,
                                edge=0.2)
        edge_sizes = system_box.lengths - packed.boundingbox.lengths
        assert np.allclose(edge_sizes, np.array([0.4]*3), atol=0.1)

        sphere = mb.fill_sphere(compound=h2o,
                                n_compounds=100,
                                sphere=[2, 2, 2, 1],
                                edge=0.2)
        assert np.allclose(sphere.boundingbox.mins, np.array([1.2]*3), atol=0.1)
        assert np.allclose(sphere.boundingbox.maxs, np.array([2.8]*3), atol=0.1)

        solvated = mb.solvate(solvent=h2o,
                              solute=methane,
                              n_solvent=100,
                              box=system_box,
                              overlap=0.2)
        edge_sizes = system_box.lengths - solvated.boundingbox.lengths
        assert np.allclose(edge_sizes, np.array([0.4]*3), atol=0.1)



        
        


