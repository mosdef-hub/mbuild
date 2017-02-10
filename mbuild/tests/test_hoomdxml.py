import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_foyer


class TestHoomdXML(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.hoomdxml')

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.hoomdxml', forcefield_name='oplsaa')

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        ethane.save(filename='ethane-box.hoomdxml', forcefield_name='oplsaa', box=box)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_rigid(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        rigid = np.zeros(ethane.n_particles)
        ethane.save(filename='ethane-box.hoomdxml',
                    forcefield_name='oplsaa',
                    box=box,
                    rigid_bodies=rigid)
