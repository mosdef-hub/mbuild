import os
import warnings

import numpy as np
import pytest

import mbuild as mb
from mbuild import Box
from mbuild.exceptions import MBuildError
from mbuild.tests.base_test import BaseTest


class TestPacking(BaseTest):
    def test_fill_box(self, h2o):
        filled = mb.fill_box(h2o, n_compounds=50, box=Box([2, 2, 2]))
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2
        assert np.array_equal(filled.box.lengths, [2, 2, 2])
        assert np.array_equal(filled.box.angles, (90, 90, 90))

    def test_fill_box_density_box(self, h2o):
        filled = mb.fill_box(h2o, n_compounds=100, density=100)
        assert np.all(
            np.isclose(filled.box.lengths, np.ones(3) * 3.104281669169261)
        )

    def test_fill_box_aspect_ratio(self, h2o):
        filled = mb.fill_box(
            h2o, n_compounds=1000, density=1000, aspect_ratio=[1, 2, 1]
        )
        assert np.isclose(filled.box.lengths[0] / filled.box.lengths[1], 0.5)
        assert np.isclose(filled.box.lengths[1] / filled.box.lengths[2], 2)

    def test_fill_box_density_n_compounds(self, h2o):
        filled = mb.fill_box(
            h2o, density=100, box=Box([3.1042931, 3.1042931, 3.1042931])
        )
        assert filled.n_particles == 300

    def test_fill_box_compound_ratio(self, h2o, ethane):
        filled = mb.fill_box(
            compound=[h2o, ethane],
            density=800,
            compound_ratio=[2, 1],
            box=Box([2, 2, 2]),
        )
        n_ethane = len([c for c in filled.children if c.name == "Ethane"])
        n_water = len([c for c in filled.children if c.name == "H2O"])
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
        filled = mb.fill_sphere(
            compound=[h2o, ethane],
            sphere=[3, 3, 3, 1.5],
            density=800,
            compound_ratio=[2, 1],
        )
        n_ethane = len([c for c in filled.children if c.name == "Ethane"])
        n_water = len([c for c in filled.children if c.name == "H2O"])
        assert n_water / n_ethane == 2

    def test_fill_sphere_bad_args(self, h2o, ethane):
        with pytest.raises(ValueError):
            mb.fill_sphere(compound=h2o, sphere=[4, 4, 4, 1])
        with pytest.raises(ValueError):
            mb.fill_sphere(
                compound=h2o, n_compounds=100, density=100, sphere=[4, 4, 4, 1]
            )
        with pytest.raises(TypeError):
            mb.fill_sphere(compound=h2o, density=1000, sphere="yes")
        with pytest.raises(ValueError):
            mb.fill_sphere(
                compound=[h2o, ethane], n_compounds=1000, sphere=[1, 1, 1, 4]
            )
        with pytest.raises(ValueError):
            mb.fill_sphere(
                compound=h2o, n_compounds=[10, 10], sphere=[1, 1, 1, 4]
            )
        with pytest.raises(ValueError):
            mb.fill_sphere(compound=h2o, n_compounds=100, sphere=[1, 1, 1, 4])

    def test_fill_region(self, h2o):
        filled = mb.fill_region(
            h2o,
            n_compounds=50,
            region=Box(lengths=[2, 3, 3], angles=[90.0, 90.0, 90.0]),
            bounds=[[3, 2, 2, 5, 5, 5]],
        )
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2
        assert np.min(filled.xyz[:, 0]) >= 3
        assert np.min(filled.xyz[:, 1]) >= 2
        assert np.min(filled.xyz[:, 2]) >= 2
        assert np.max(filled.xyz[:, 0]) <= 5
        assert np.max(filled.xyz[:, 1]) <= 5
        assert np.max(filled.xyz[:, 2]) <= 5

    def test_fill_region_box(self, h2o):
        mybox = Box(lengths=[4, 4, 4], angles=[90.0, 90.0, 90.0])
        filled = mb.fill_region(
            h2o, n_compounds=50, region=mybox, bounds=[[0, 0, 0, 4, 4, 4]]
        )
        assert filled.n_particles == 50 * 3
        assert filled.n_bonds == 50 * 2
        assert np.min(filled.xyz[:, 0]) >= 0
        assert np.max(filled.xyz[:, 2]) <= 4

    def test_fill_region_multiple(self, ethane, h2o):
        box1 = mb.Box(lengths=[2, 2, 2], angles=[90.0, 90.0, 90.0])
        box2 = mb.Box(lengths=[2, 2, 2], angles=[90.0, 90.0, 90.0])
        filled = mb.fill_region(
            compound=[ethane, h2o],
            n_compounds=[2, 2],
            region=[box1, box2],
            bounds=[[2, 2, 2, 4, 4, 4], [4, 2, 2, 6, 4, 4]],
        )
        assert filled.n_particles == 2 * 8 + 2 * 3
        assert filled.n_bonds == 2 * 7 + 2 * 2
        assert np.max(filled.xyz[:16, 0]) < 4
        assert np.min(filled.xyz[16:, 0]) > 4

    def test_fill_region_incorrect_type(self, ethane):
        box1 = {"a": 1}
        with pytest.raises(ValueError, match=r"expected a list of type:"):
            mb.fill_region(
                compound=[ethane], n_compounds=[2], region=box1, bounds=None
            )

    def test_fill_region_bounds_not_list(self, ethane):
        box1 = Box(lengths=[2, 2, 2], angles=[90.0, 90.0, 90.0])
        with pytest.raises(TypeError):
            mb.fill_region(
                compound=[ethane], n_compounds=[2], region=box1, bounds=box1
            )

    def test_fill_region_incorrect_bounds_amount(self, ethane, h2o):
        box1 = Box(lengths=[2, 2, 2], angles=[90.0, 90.0, 90.0])
        with pytest.raises(ValueError):
            mb.fill_region(
                compound=[ethane, h2o],
                n_compounds=[2, 2],
                region=box1,
                bounds=[box1],
            )

    def test_fill_region_incorrect_bounds_types(self, ethane, h2o):
        box1 = Box(lengths=[2, 2, 2], angles=[90.0, 90.0, 90.0])
        with pytest.raises(ValueError):
            mb.fill_region(
                compound=[ethane, h2o],
                n_compounds=[2, 2],
                region=box1,
                bounds=[1.0, 1.0, 1.0],
            )

    def test_box_no_bound(self, ethane):
        box1 = Box(lengths=[2, 2, 2], angles=[90.0, 90.0, 90.0])
        mb.fill_region(
            compound=[ethane], n_compounds=[2], region=box1, bounds=None
        )

    def test_fill_region_multiple_bounds(self, ethane, h2o):
        box1 = Box.from_mins_maxs_angles(
            mins=[2, 2, 2], maxs=[4, 4, 4], angles=[90.0, 90.0, 90.0]
        )
        box2 = mb.Box.from_mins_maxs_angles(
            mins=[4, 2, 2], maxs=[6, 4, 4], angles=[90.0, 90.0, 90.0]
        )
        filled = mb.fill_region(
            compound=[ethane, h2o],
            n_compounds=[2, 2],
            region=[box1, box2],
            bounds=[[2, 2, 2, 4, 4, 4], [4, 2, 2, 6, 4, 4]],
        )
        assert filled.n_particles == 2 * 8 + 2 * 3
        assert filled.n_bonds == 2 * 7 + 2 * 2
        assert np.max(filled.xyz[:16, 0]) < 4
        assert np.min(filled.xyz[16:, 0]) > 4

    def test_fill_region_multiple_types(self, ethane, h2o):
        box1 = mb.Box.from_mins_maxs_angles(
            mins=[2, 2, 2], maxs=[4, 4, 4], angles=[90.0, 90.0, 90.0]
        )
        box2 = [4, 2, 2, 6, 4, 4]
        filled = mb.fill_region(
            compound=[ethane, h2o],
            n_compounds=[2, 2],
            region=[box1, box2],
            bounds=[[2, 2, 2, 4, 4, 4], box2],
        )
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
        ethane_pos = ethane.pos
        solvated1 = mb.solvate(
            ethane, h2o, n_solvent=n_solvent, box=[4, 4, 4], center_solute=True
        )
        solvated2 = mb.solvate(
            ethane, h2o, n_solvent=n_solvent, box=[4, 4, 4], center_solute=False
        )

        print(ethane.xyz)
        print(solvated2.children[0].xyz)

        assert (
            solvated1.n_particles == solvated2.n_particles == 8 + n_solvent * 3
        )
        assert solvated1.n_bonds == solvated2.n_bonds == 7 + n_solvent * 2
        assert not np.isclose(solvated1.children[0].pos, ethane_pos).all()
        assert np.isclose(solvated2.children[0].pos, ethane_pos).all()

        assert ethane.parent == None

    def test_solvate_multiple(self, methane, ethane, h2o):
        init_box = mb.fill_box(methane, 2, box=[4, 4, 4])
        solvated = mb.solvate(init_box, [ethane, h2o], [20, 20], box=[4, 4, 4])
        assert solvated.n_particles == 2 * 5 + 20 * 8 + 20 * 3
        assert len(solvated.children) == 41

    def test_fill_box_seed(self, ethane):
        filled = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2])
        filled_same = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2])
        filled_diff = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2], seed=2)
        assert np.array_equal(filled.xyz, filled_same.xyz)
        assert not np.array_equal(filled.xyz, filled_diff.xyz)

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
            mb.solvate(
                solute=h2o, solvent=[h2o], n_solvent=[10, 10], box=[2, 2, 2]
            )
        with pytest.raises(ValueError):
            mb.fill_region(h2o, n_compounds=[10, 10], region=[2, 2, 2, 4, 4, 4])
        with pytest.raises(ValueError):
            mb.fill_box(
                compound=[h2o, h2o],
                n_compounds=[10],
                density=1000,
                fix_orientation=[True, True, True],
            )

    def test_write_temp_file(self, h2o):
        cwd = os.getcwd()  # Must keep track of the temp dir that pytest creates
        filled = mb.fill_box(
            h2o, n_compounds=10, box=Box([4, 4, 4]), temp_file="temp_file1.pdb"
        )
        region = mb.fill_region(
            h2o,
            10,
            [[2, 2, 2, 4, 4, 4]],
            temp_file="temp_file2.pdb",
            bounds=[[2, 2, 2, 4, 4, 4]],
        )
        solvated = mb.solvate(
            filled, h2o, 10, box=[4, 4, 4], temp_file="temp_file3.pdb"
        )
        assert os.path.isfile(os.path.join(cwd, "temp_file1.pdb"))
        assert os.path.isfile(os.path.join(cwd, "temp_file2.pdb"))
        assert os.path.isfile(os.path.join(cwd, "temp_file3.pdb"))

    def test_packmol_error(self, h2o):
        with pytest.raises(MBuildError, match=r"co\-linear"):
            mb.fill_box(h2o, n_compounds=10, box=[0, 0, 0])

    def test_packmol_warning(self, h2o):
        with pytest.warns(UserWarning):
            mb.fill_box(h2o, n_compounds=10, box=[1, 1, 1], overlap=10)

    def test_packmol_args(self, h2o):
        try:
            mb.fill_box(
                h2o,
                n_compounds=2,
                box=[0.1, 0.1, 0.1],
                packmol_args={"maxit": 10, "movebadrandom": "", "nloop": 100},
            )
        except RuntimeError:
            with open("log.txt", "r") as logfile:
                assert "(movebadrandom)" in logfile.read()
                logfile.seek(0)
                assert (
                    "User defined GENCAN number of iterations:           10"
                    in logfile.read()
                )

        try:
            mb.fill_region(
                h2o,
                10,
                [[2, 2, 2, 4, 4, 4]],
                bounds=[[2, 2, 2, 4, 4, 4]],
                packmol_args={"maxit": 10, "movebadrandom": "", "nloop": 100},
            )
        except RuntimeError:
            with open("log.txt", "r") as logfile:
                assert "(movebadrandom)" in logfile.read()
                logfile.seek(0)
                assert (
                    "User defined GENCAN number of iterations:           10"
                    in logfile.read()
                )

        try:
            mb.fill_sphere(
                h2o,
                n_compounds=2,
                sphere=[3, 3, 3, 1.5],
                packmol_args={"maxit": 10, "movebadrandom": "", "nloop": 100},
            )
        except RuntimeError:
            with open("log.txt", "r") as logfile:
                assert "(movebadrandom)" in logfile.read()
                logfile.seek(0)
                assert (
                    "User defined GENCAN number of iterations:           10"
                    in logfile.read()
                )
        try:
            mb.solvate(
                solute=h2o,
                solvent=[h2o],
                n_solvent=[10],
                box=[2, 2, 2],
                packmol_args={"maxit": 15, "movebadrandom": "", "nloop": 100},
            )
        except RuntimeError:
            with open("log.txt", "r") as logfile:
                assert "(movebadrandom)" in logfile.read()
                logfile.seek(0)
                assert (
                    "User defined GENCAN number of iterations:           15"
                    in logfile.read()
                )

    def test_packmol_args_bad(self, h2o):
        with pytest.raises(ValueError):
            mb.fill_box(
                h2o,
                n_compounds=2,
                box=[10, 10, 10],
                packmol_args={"bad_arg": 10},
            )

    @pytest.mark.parametrize(
            "args",
            [
                dict(maxit=500),
                dict(nloop=1000),
                dict(movebadrandom=""),
                dict(fbins=1.2),
                dict(discale=1.5),
                dict(movefrac=0.05),
                dict(avoid_overlap=""),
                dict(precision=0.02),
                dict(use_short_tol="", short_tol_dist=0.5),
                dict(use_short_tol="", short_tol_scale=1.2)
            ]
    )
    def test_packmol_args_allowed(self, args):
            mb.fill_box(
                mb.load("C", smiles=True),
                n_compounds=10,
                box=[10, 10, 10],
                packmol_args=args,
            )

    @pytest.mark.parametrize(
            "args",
            [
                dict(tolerance=0.2),
                dict(seed=42),
                dict(sidemax=2.0),
            ]
    )
    def test_packmol_args_default(self, args):
        with pytest.warns():
            mb.fill_box(
                mb.load("C", smiles=True),
                n_compounds=10,
                box=[10, 10, 10],
                packmol_args=args,
            )


    def test_rotate(self, h2o):
        filled = mb.fill_box(h2o, 2, box=[1, 1, 1], fix_orientation=True)
        w0 = filled.xyz[:3]
        w1 = filled.xyz[3:]
        # Translate w0 and w1 to COM
        w0 -= w0.sum(0) / len(w0)
        w1 -= w1.sum(0) / len(w1)
        assert np.isclose(w0, w1).all()

    def test_no_rotate(self, h2o):
        filled = mb.fill_box(
            [h2o, h2o], [1, 1], box=[1, 1, 1], fix_orientation=[False, True]
        )
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
        # With default sidemax
        box_of_methane = mb.fill_box(
            ch4, box=[1000, 1000, 1000], n_compounds=500
        )
        sphere_of_methane = mb.fill_sphere(
            ch4, sphere=[1000, 1000, 1000, 1000], n_compounds=500
        )
        assert all(
            np.asarray(box_of_methane.get_boundingbox().lengths)
            < [110, 110, 110]
        )
        assert all(
            np.asarray(sphere_of_methane.get_boundingbox().lengths)
            < [210, 210, 210]
        )

        # With adjusted sidemax
        big_box_of_methane = mb.fill_box(
            ch4, box=[1000, 1000, 1000], n_compounds=500, sidemax=1000.0
        )
        big_sphere_of_methane = mb.fill_sphere(
            ch4,
            sphere=[1000, 1000, 1000, 1000],
            n_compounds=500,
            sidemax=2000.0,
        )

        assert all(
            np.asarray(big_box_of_methane.get_boundingbox().lengths)
            > [900, 900, 900]
        )
        assert all(
            np.asarray(big_sphere_of_methane.get_boundingbox().lengths)
            > [1800, 1800, 1800]
        )

    def test_box_edge(self, h2o, methane):
        system_box = mb.Box(lengths=(1.8, 1.8, 1.8))
        packed = mb.fill_box(
            compound=h2o, n_compounds=100, box=system_box, edge=0.2
        )
        edge_sizes = np.subtract(
            system_box.lengths, packed.get_boundingbox().lengths
        )
        assert np.allclose(edge_sizes, np.array([0.4] * 3), atol=0.1)

        region = mb.fill_region(
            compound=h2o,
            n_compounds=100,
            region=system_box,
            edge=0.2,
            bounds=[system_box],
        )
        edge_sizes = np.subtract(
            system_box.lengths, packed.get_boundingbox().lengths
        )
        assert np.allclose(edge_sizes, np.array([0.4] * 3), atol=0.1)

        edge = 0.2
        bounds = [2, 2, 2, 1]
        sphere = mb.fill_sphere(
            compound=h2o, n_compounds=120, sphere=bounds, edge=edge
        )
        target_diameter = (bounds[3] - edge) * 2
        assert np.allclose(
            sphere.maxs - sphere.mins, np.array([target_diameter] * 3), atol=0.1
        )

        solvated = mb.solvate(
            solvent=h2o,
            solute=methane,
            n_solvent=100,
            box=system_box,
            overlap=0.2,
        )
        edge_sizes = np.subtract(
            system_box.lengths, solvated.get_boundingbox().lengths
        )
        assert np.allclose(edge_sizes, np.array([0.4] * 3), atol=0.1)

    def test_validate_mass(self, methane):
        bead = mb.Compound(name="A", mass=0.0)
        with pytest.raises(MBuildError):
            mb.fill_box(compound=bead, n_compounds=10, density=1)

        with pytest.raises(MBuildError):
            mb.fill_box(compound=bead, density=1, box=[0.5, 0.5, 0.5])

        with pytest.raises(MBuildError):
            mb.fill_sphere(compound=bead, sphere=[20, 20, 20, 20], density=1)

        beadA = mb.Compound(name="A", mass=0.0)
        beadB = mb.Compound(name="B", mass=1.0, pos=[0.5, 0.5, 0.5])
        beads = mb.Compound(subcompounds=[beadA, beadB])
        with warnings.catch_warnings(record=True) as w:
            mb.packing._validate_mass(compound=[beadA, beadB], n_compounds=None)
            assert w

        with warnings.catch_warnings(record=True) as w:
            mb.packing._validate_mass(compound=[beads], n_compounds=None)
            assert w

        with warnings.catch_warnings(record=True) as w:
            mb.packing._validate_mass(compound=[beads], n_compounds=[5])
            assert w
