import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestXYZ(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.xyz')
        ethane.save(filename='ethane.pdb')
        ethane_in = mb.load('ethane.xyz', top='ethane.pdb')
        assert len(ethane_in.children) == 8
        assert set([child.name for child in ethane_in.children]) == {'C', 'H'}
