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

    def test_dtype(self):
        box = mb.Box(mins=np.zeros(3), maxs=np.ones(3))
        assert box.lengths.dtype == np.float64
        assert box.mins.dtype == np.float64
        assert box.maxs.dtype == np.float64

    def test_mins_setter(self):
        box = mb.Box(mins=np.zeros(3), maxs=2 * np.ones(3))
        box.mins = np.ones(3)
        assert (box.mins == np.ones(3)).all()
        assert (box.maxs - box.mins == np.ones(3)).all()
        assert (box.lengths == np.ones(3)).all()

    def test_maxs_setter(self):
        box = mb.Box(mins=np.zeros(3), maxs=np.ones(3))
        box.maxs = 2 * np.ones(3)
        assert (box.maxs == 2 * np.ones(3)).all()
        assert (box.maxs - box.mins == 2 * np.ones(3)).all()
        assert (box.lengths == 2 * np.ones(3)).all()

    def test_lengths_setter(self):
        box = mb.Box(mins=np.zeros(3), maxs=np.ones(3))
        box.lengths = 2 * np.ones(3)
        assert (box.lengths == 2 * np.ones(3)).all()
        assert (box.maxs - box.mins == 2 * np.ones(3)).all()
