import numpy as np
import pytest

import mbuild as mb
import mbuild.lib.molecules.water as water_models
from mbuild.exceptions import MBuildError
from mbuild.lib.recipes.water_box import Water3SiteBox
from mbuild.tests.base_test import BaseTest


class TestWaterBox(BaseTest):
    def test_replicate(self):
        temp_box = mb.Box([1.0, 1.0, 1.0])
        wb = Water3SiteBox(box=temp_box)
        assert wb.n_particles == 72

        temp_box = mb.Box([1.0, 1.0, 1.0])
        wb = Water3SiteBox(box=temp_box, edge=0.15)
        assert wb.n_particles == 60

        temp_box = mb.Box([2.0, 2.0, 2.0])
        wb = Water3SiteBox(box=temp_box)
        assert wb.n_particles == 678

        wb = Water3SiteBox(box=temp_box, edge=[0.2, 0.2, 0.2])
        assert wb.n_particles == 567

        wb = Water3SiteBox(box=temp_box, edge=0.2)
        assert wb.n_particles == 567

        wb = Water3SiteBox(box=temp_box, edge=[0.2, 0.2, 0.5])
        assert wb.n_particles == 474

        temp_box = mb.Box([8.0, 1.0, 1.0])
        wb = Water3SiteBox(box=temp_box)
        assert wb.n_particles == 657
        assert wb.box == temp_box

        temp_box = mb.Box([1.0, 8.0, 1.0])
        wb = Water3SiteBox(box=temp_box)
        assert wb.n_particles == 606
        assert wb.box == temp_box

        temp_box = mb.Box([1.0, 1.0, 8.0])
        wb = Water3SiteBox(box=temp_box)
        assert wb.n_particles == 603
        assert wb.box == temp_box

        wb = Water3SiteBox(box=mb.Box([2.0, 5.5, 13.2]))
        assert wb.n_particles == 13227

    def test_mask(self, ethane):
        box_temp = mb.Box([2.0, 2.0, 1.0])

        ethane_system = mb.fill_box(
            ethane,
            n_compounds=50,
            overlap=0.22,
            edge=0.10,
            sidemax=box_temp.Lz + 1,
            box=box_temp,
            seed=1234,
        )

        wb = Water3SiteBox(box=mb.Box([2.0, 2.0, 2.0]), mask=ethane_system)

        assert wb.n_particles == 285

        wb = Water3SiteBox(
            box=mb.Box([2.0, 2.0, 2.0]), mask=ethane_system, radii_scaling=0.9
        )

        assert wb.n_particles == 294

        wb = Water3SiteBox(
            box=mb.Box([2.0, 2.0, 2.0]), mask=[ethane_system], radii_scaling=0.9
        )

        assert wb.n_particles == 294

        wb = Water3SiteBox(
            box=mb.Box([2.0, 2.0, 2.0]), mask=[ethane_system], radii_scaling=1.1
        )

        assert wb.n_particles == 282

        # test using the default cutoff for ethane
        for part in ethane.particles():
            part.element = None

        ethane_system = mb.fill_box(
            ethane,
            n_compounds=50,
            overlap=0.22,
            edge=0.10,
            sidemax=box_temp.Lz + 1,
            box=box_temp,
            seed=1234,
        )

        wb = Water3SiteBox(
            box=mb.Box([2.0, 2.0, 2.0]), mask=ethane_system, radii_overlap=0.1
        )

        assert wb.n_particles == 294

        wb = Water3SiteBox(
            box=mb.Box([2.0, 2.0, 2.0]), mask=ethane_system, radii_overlap=0.2
        )

        assert wb.n_particles == 249

        wb = Water3SiteBox(
            box=mb.Box([2.0, 2.0, 2.0]),
            mask=ethane_system,
            radii_dict={"C": 0.1, "H": 0.05},
        )

        assert wb.n_particles == 321

        wb = Water3SiteBox(
            box=mb.Box([2.0, 2.0, 2.0]),
            mask=ethane_system,
            radii_dict={"C": 0.15, "H": 0.125},
        )

        assert wb.n_particles == 285

    def test_model(self):
        wb = Water3SiteBox(
            box=mb.Box([2.0, 2.0, 2.0]), model=water_models.WaterSPC()
        )
        assert wb.children[0].name == "WaterSPC"

        parts = [p for p in wb.children[0].particles()]
        oh_dist1 = np.linalg.norm(parts[0].pos - parts[1].pos)
        oh_dist2 = np.linalg.norm(parts[0].pos - parts[2].pos)
        assert np.isclose(oh_dist1, 0.1)
        assert np.isclose(oh_dist2, 0.1)

        v1 = np.array(parts[1].pos) - np.array(parts[0].pos)
        v2 = np.array(parts[2].pos) - np.array(parts[0].pos)
        v1_u = v1 / np.linalg.norm(v1)

        v2_u = v2 / np.linalg.norm(v2)
        angle = (
            np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) * 180.0 / np.pi
        )

        assert np.isclose(angle, 109.47)

    def test_incorrect_edge(self):
        with pytest.raises(ValueError):
            Water3SiteBox(box=mb.Box([2.0, 2.0, 2.0]), edge=[1,2,3,4])

    def test_incorrect_edge2(self):
        with pytest.raises(ValueError):
            Water3SiteBox(box=mb.Box([2.0, 2.0, 2.0]), edge=[1, 2])

    def test_incorrect_mask(self):
        with pytest.raises(MBuildError):
            Water3SiteBox(box=mb.Box([2.0, 2.0, 2.0]), mask=[1,2,3,4])

    def test_incorrect_mask2(self):
        with pytest.raises(MBuildError):
            Water3SiteBox(box=mb.Box([2.0, 2.0, 2.0]), mask=1)

    def test_bad_model(self):
        bad_model = mb.Compound()
        oxygen = mb.Compound(name="O", element="O")
        hydrogen = mb.Compound(name="H", element="H")
        bad_model.add(
            [mb.clone(hydrogen), mb.clone(oxygen), mb.clone(hydrogen)]
        )

        with pytest.raises(MBuildError):
            Water3SiteBox(box=mb.Box([2.0, 2.0, 2.0]), model=bad_model)

    def test_bad_model2(self):
        bad_model = mb.Compound()
        oxygen = mb.Compound(name="O", element="O")
        hydrogen = mb.Compound(name="H", element="H")
        bad_model.add(
            [mb.clone(hydrogen), mb.clone(oxygen), mb.clone(hydrogen)]
        )

        with pytest.raises(MBuildError):
            Water3SiteBox(box=mb.Box([2.0, 2.0, 2.0]), model=[bad_model])

    def test_bad_dict(self):
        bad_model = mb.Compound()
        oxygen = mb.Compound(name="O", element="O")
        hydrogen = mb.Compound(name="H", element="H")
        bad_model.add(
            [mb.clone(hydrogen), mb.clone(oxygen), mb.clone(hydrogen)]
        )

        with pytest.raises(ValueError):
            Water3SiteBox(box=mb.Box([2.0, 2.0, 2.0]), mask=oxygen, radii_dict=[1.0,2.0])
