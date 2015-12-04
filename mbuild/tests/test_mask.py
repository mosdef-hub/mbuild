import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestMask(BaseTest):

    @pytest.mark.skipif(True, reason='Needs implementing!')
    def test_apply_pattern(self):
        pass

    def test_random_2d(self):
        pattern = mb.Random2DPattern(100)
        assert len(pattern) == 100

    def test_random_3d(self):
        pattern = mb.Random3DPattern(100)
        assert len(pattern) == 100

    def test_grid_2d(self):
        pattern = mb.Grid2DPattern(10, 5)
        assert len(pattern) == 50

    def test_grid_3d(self):
        pattern = mb.Grid3DPattern(10, 5, 2)
        assert len(pattern) == 100

    def test_sphere(self):
        pattern = mb.SpherePattern(100)
        assert len(pattern) == 100

    def test_disk(self):
        pattern = mb.DiskPattern(100)
        assert len(pattern) == 100
