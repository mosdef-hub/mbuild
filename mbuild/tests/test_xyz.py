import numpy as np
import pytest

import mbuild as mb
from mbuild.utils.io import get_fn
from mbuild.tests.base_test import BaseTest
from mbuild.exceptions import MBuildError


class TestXYZ(BaseTest):
    def test_load_no_top(self, ethane):
        ethane.save(filename='ethane.xyz')
        ethane_in = mb.load('ethane.xyz')
        assert len(ethane_in.children) == 8
        assert ethane_in.n_bonds == 0
        assert set([child.name for child in ethane_in.children]) == {'C', 'H'}

    def test_wrong_n_atoms(self):
        with pytest.raises(MBuildError):
            mb.load(get_fn('too_few_atoms.xyz'))
        with pytest.raises(MBuildError):
            mb.load(get_fn('too_many_atoms.xyz'))

    def test_save(self, ethane):
        ethane.save(filename='ethane.xyz')
        ethane_in = mb.load('ethane.xyz')
        assert len(ethane_in.children) == 8
        assert set([child.name for child in ethane_in.children]) == {'C', 'H'}

    def test_coordinates(self, ethane):
        ethane.save(filename='ethane.xyz')
        ethane_in = mb.load('ethane.xyz')
        assert np.allclose(ethane.xyz, ethane_in.xyz, atol=1e-3)
