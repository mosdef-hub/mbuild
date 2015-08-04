import pytest
import sys

import mbuild as mb
from mbuild.tests.base_test import BaseTest

class TestPacking(BaseTest):

    def test_fill_box(self, ethane):
        filled = mb.fill_box(ethane, n_compounds=20, box=[2, 2, 2])
        assert filled.n_atoms == 20 * 8
        assert filled.n_bonds == 20 * 7

    def test_solvate(self, ethane, h2o):
        solvated = mb.solvate(ethane, h2o, n_solvent=100, box=[4, 4, 4])
        assert solvated.n_atoms == 8 + 100 * 3
        assert solvated.n_bonds == 7 + 100 * 2
