import ele
import numpy as np
import pytest

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.formats.xyz import write_xyz
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn


class TestXYZ(BaseTest):
    def test_load_no_top(self, ethane):
        ethane.save(filename="ethane.xyz")
        ethane_in = mb.load("ethane.xyz")
        assert len(ethane_in.children) == 1
        assert ethane_in.n_bonds == 0
        assert set([p.name for p in ethane_in.particles()]) == {"C", "H"}
        assert set([p.element for p in ethane_in.particles()]) == {
            ele.Elements.C,
            ele.Elements.H,
        }

    def test_wrong_n_atoms(self):
        with pytest.raises(MBuildError):
            mb.load(get_fn("too_few_atoms.xyz"))
        with pytest.raises(MBuildError):
            mb.load(get_fn("too_many_atoms.xyz"))

    def test_bad_input(self, ethane):
        with pytest.raises(ValueError):
            assert isinstance(ethane, mb.Compound)
            write_xyz(ethane, "compound.xyz")

    def test_save(self, ethane):
        ethane.save(filename="ethane.xyz")
        ethane_in = mb.load("ethane.xyz")
        assert len(ethane_in.children) == 8
        assert ethane_in.n_bonds == 0
        assert set([child.name for child in ethane_in.children]) == {"C", "H"}

    def test_coordinates(self, ethane):
        ethane.save(filename="ethane.xyz")
        ethane_in = mb.load("ethane.xyz")
        assert np.allclose(ethane.xyz, ethane_in.xyz)

    def test_non_resolved_elements(self):
        tip3p_water = mb.load(get_fn("tip3p_water.xyz"))
        assert tip3p_water[0].element is None
        assert tip3p_water[0].name == "opls_111"
        assert tip3p_water[1].element is None
        assert tip3p_water[1].name == "opls_112"
        assert tip3p_water[2].element is None
        assert tip3p_water[2].name == "opls_112"

    def test_write_atomtypes(self):
        tip3p_water = mb.load(get_fn("tip3p_water.xyz"))
        tip3p_water.save(filename="test.xyz", write_atomnames=True)
        lines = []
        with open("test.xyz") as f:
            for line in f:
                lines.append(line.split())
        assert lines[2][0] == "opls_111"
        assert lines[3][0] == "opls_112"
        assert lines[4][0] == "opls_112"
        tip3p_water_in = mb.load("test.xyz")
        assert tip3p_water_in[0].element is None
        assert tip3p_water_in[0].name == "opls_111"
        assert tip3p_water_in[1].element is None
        assert tip3p_water_in[1].name == "opls_112"
        assert tip3p_water_in[2].element is None
        assert tip3p_water_in[2].name == "opls_112"
