import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestMask(BaseTest):

    @pytest.mark.skipif(True, reason='Needs implementing!')
    def test_apply_mask(self):
        pass

    def test_random_2d(self):
        mask = mb.random_mask_2d(100)
        assert len(mask) == 100

    def test_random_3d(self):
        mask = mb.random_mask_3d(100)
        assert len(mask) == 100

    def test_grid_2d(self):
        mask = mb.grid_mask_2d(10, 5)
        assert len(mask) == 50

    def test_grid_3d(self):
        mask = mb.grid_mask_3d(10, 5, 2)
        assert len(mask) == 100

    def test_sphere(self):
        mask = mb.sphere_mask(100)
        assert len(mask) == 100

    def test_disk(self):
        mask = mb.disk_mask(100)
        assert len(mask) == 100
