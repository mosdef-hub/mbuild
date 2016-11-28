import pytest
import numpy as np
import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestPattern(BaseTest):

    def test_apply_to_compound(self, betacristobalite, alkyl, ch3):
        pattern = mb.Random2DPattern(90)
        chains, backfills = pattern.apply_to_compound(
            guest=alkyl, host=betacristobalite, backfill=ch3)
        assert len(chains) == 90
        assert len(backfills) == 10

        with pytest.raises(AssertionError):
            pattern = mb.Random2DPattern(101)
            chains, backfills = pattern.apply_to_compound(
                guest=alkyl, host=betacristobalite, backfill=ch3)

    def test_random_2d(self):
        pattern = mb.Random2DPattern(100)
        assert len(pattern) == 100

    def test_random_2d_seed(self):
        pattern_a = mb.Random2DPattern(100, seed=12345)
        pattern_b = mb.Random2DPattern(100, seed=12345)
        assert np.array_equal(pattern_a, pattern_b)

    def test_random_3d(self):
        pattern = mb.Random3DPattern(100)
        assert len(pattern) == 100

    def test_random_3d_seed(self):
        pattern_a = mb.Random3DPattern(100, seed=12345)
        pattern_b = mb.Random3DPattern(100, seed=12345)
        assert np.array_equal(pattern_a, pattern_b)

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
