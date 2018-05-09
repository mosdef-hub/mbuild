import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestXYZ(BaseTest):
    def test_load_no_top(self, ethane):
        ethane.save(filename='ethane.xyz')
        ethane.save(filename='ethane.pdb')
        ethane_in = mb.load('ethane.xyz')
        assert len(ethane_in.children) == 8
        assert set([child.name for child in ethane_in.children]) == {'C', 'H'}

    def test_load_with_top(self, ethane):
        ethane.save(filename='ethane.xyz')
        ethane.save(filename='ethane.pdb')
        ethane_in = mb.load('ethane.xyz', top='ethane.pdb')
        assert len(ethane_in.children) == 8
        ethane_in.n_bonds == 7
        assert set([child.name for child in ethane_in.children]) == {'C', 'H'}
