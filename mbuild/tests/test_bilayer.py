import numpy as np
import math
import pytest
from mbuild.tests.base_test import BaseTest
from mbuild.recipes.bilayer.bilayer import Bilayer

from mbuild.exceptions import MBuildError


class TestBilayer(BaseTest):
    """
    Unit Tests for Bilayer class functionality.
    """
    def test_size(self, bilayer_default, dspcua):
        assert len(bilayer_default['lipid_bilayer']['top_leaflet'].children) == 4
        assert len(bilayer_default['lipid_bilayer']['bottom_leaflet'].children) == 4
        assert bilayer_default['lipid_bilayer']['top_leaflet'].n_particles == dspcua.n_particles * 4
        assert bilayer_default['lipid_bilayer']['bottom_leaflet'].n_particles == dspcua.n_particles * 4
        assert bilayer_default['solvent'].children[0].n_particles / 3 == \
               20 * bilayer_default.n_lipids_x * bilayer_default.n_lipids_y * 2

    def test_composition(self, ternary_lipid_mix):
        bilayer = Bilayer(ternary_lipid_mix, n_lipids_x=2, n_lipids_y=2, solvent_per_lipid=0)
        assert bilayer.number_of_each_lipid_per_leaflet == [2, 1, 1]
        bilayer = Bilayer(ternary_lipid_mix, n_lipids_x=4, n_lipids_y=2, solvent_per_lipid=0)
        assert bilayer.number_of_each_lipid_per_leaflet == [4, 2, 2]

    def test_leaflet_spacing(self, bilayer_default):
        top_z_min = np.amin(bilayer_default['lipid_bilayer']['top_leaflet'].xyz, axis=0)[2]
        bottom_z_max = np.amax(bilayer_default['lipid_bilayer']['bottom_leaflet'].xyz, axis=0)[2]
        space = top_z_min - bottom_z_max
        assert math.isclose(space, bilayer_default.z_spacing)

    def test_bad_dimensions(self, binary_lipid_mix, ternary_lipid_mix):
        with pytest.raises(ValueError):
            Bilayer(binary_lipid_mix, n_lipids_x=0, n_lipids_y=1)
        with pytest.raises(TypeError):
            Bilayer(binary_lipid_mix, n_lipids_x=4.5, n_lipids_y=5)
        with pytest.raises(ValueError):
            Bilayer(ternary_lipid_mix, n_lipids_x=2, n_lipids_y=1)

    def test_bad_lipids(self, dspcua, ffaua, alcua):
        with pytest.raises(TypeError):
            Bilayer((dspcua, ffaua, alcua, dspcua))
        with pytest.raises(ValueError):
            Bilayer((dspcua, 1.0, 0), n_lipids_x=2, n_lipids_y=2)
        with pytest.raises(ValueError):
            lipids = [(ffaua, 0.5, 0.0, 0), (alcua, 0.5, 17)]
            Bilayer(lipids, n_lipids_x=2, n_lipids_y=2)
        with pytest.raises(ValueError):
            frac = [(ffaua, 0.7, 0.0, 0), (alcua, 0.5, 17)]
            Bilayer(frac, n_lipids_x=2, n_lipids_y=2)
        with pytest.raises(ValueError):
            Bilayer(lipids, n_lipids_x=1, n_lipids_y=1)

    def test_bad_composition(self, dspcua, ffaua, alcua):
        with pytest.raises(MBuildError):
            lipids = [(dspcua, 0.33, 0.0, 0), (alcua, 0.33, 0.0, 17),
                      (ffaua, 0.34, 0.0, 17)]
            Bilayer(lipids, n_lipids_x=2, n_lipids_y=2)

    def test_bad_tilt_angle(self, binary_lipid_mix):
        with pytest.raises(TypeError):
            Bilayer(binary_lipid_mix, tilt_angle="35")

    def test_bad_random(self, binary_lipid_mix):
        with pytest.raises(TypeError):
            Bilayer(binary_lipid_mix, max_tail_randomization="24")
        with pytest.raises(ValueError):
            Bilayer(binary_lipid_mix, max_tail_randomization=150)

    def test_bad_solvent(self, binary_lipid_mix):
        with pytest.raises(TypeError):
            Bilayer(binary_lipid_mix, solvent=5)
        with pytest.raises(ValueError):
            Bilayer(binary_lipid_mix, solvent_per_lipid=10.5)

    def test_bad_itp(self, binary_lipid_mix):
        with pytest.raises(TypeError):
            Bilayer(binary_lipid_mix, itp_path=24)
        with pytest.raises(IOError):
            Bilayer(binary_lipid_mix, itp_path="bruh")

    def test_bad_files(self, binary_lipid_mix):
        with pytest.raises(TypeError):
            Bilayer(binary_lipid_mix, make_files="true")
        with pytest.raises(TypeError):
            Bilayer(binary_lipid_mix, filename=67)