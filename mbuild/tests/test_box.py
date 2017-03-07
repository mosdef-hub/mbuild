import pytest
import numpy as np
import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestBox(BaseTest):

    def test_init_lengths(self):
        box = mb.Box(lengths=np.ones(3))
        assert np.array_equal(box.lengths, np.ones(3))
        assert np.array_equal(box.mins, np.zeros(3))
        assert np.array_equal(box.maxs, np.ones(3))

    def test_init_bounds(self):
        box = mb.Box(mins=np.zeros(3), maxs=np.ones(3))
        assert np.array_equal(box.lengths, np.ones(3))
        assert np.array_equal(box.mins, np.zeros(3))
        assert np.array_equal(box.maxs, np.ones(3))

    def test_scale(self):
        box = mb.Box(lengths=np.ones(3))
        scaling_factors = np.array([3, 4, 5])
        box.scale(scaling_factors)
        assert np.array_equal(box.lengths, scaling_factors)
        assert np.array_equal(box.mins, (np.ones(3) / 2) - (scaling_factors / 2))
        assert np.array_equal(box.maxs, (scaling_factors / 2) + (np.ones(3) / 2))

    def test_center(self):
        box = mb.Box(lengths=np.ones(3))
        box.center()
        assert np.array_equal(box.lengths, np.ones(3))
        assert np.array_equal(box.mins, np.ones(3) * -0.5)
        assert np.array_equal(box.maxs, np.ones(3) * 0.5)
