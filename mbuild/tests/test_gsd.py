import mbuild as mb
import numpy as np
import pytest
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_foyer


class TestGSD(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.gsd')

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.gsd',forcefield='opls')

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0,2.0,2.0]))
        ethane.save(filename='ethane-box.gsd',forcefield='opls',box=box)
