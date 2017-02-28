import mbuild as mb
import numpy as np
import pytest
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_foyer, has_gsd


class TestGSD(BaseTest):

    @pytest.mark.skipif(not has_gsd, reason="GSD package not installed")
    def test_save(self, ethane):
        ethane.save(filename='ethane.gsd')

    @pytest.mark.skipif(not has_gsd, reason="GSD package not installed")
    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.gsd', forcefield_name='oplsaa')

    @pytest.mark.skipif(not has_gsd, reason="GSD package not installed")
    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0,2.0,2.0]))
        ethane.save(filename='ethane-box.gsd', forcefield_name='oplsaa',box=box)
