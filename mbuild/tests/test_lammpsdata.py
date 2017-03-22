import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestLammpsData(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.lammps')

    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.lammps', forcefield_name='oplsaa')

    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        ethane.save(filename='ethane-box.lammps', forcefield_name='oplsaa', box=box)
