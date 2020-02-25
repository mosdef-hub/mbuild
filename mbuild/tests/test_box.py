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

    def test_bad_args(self):
        with pytest.raises(ValueError):
            mb.Box(maxs=[4, 4, 4])

    def test_init_angles(self):
        box = mb.Box(mins=np.zeros(3), maxs=np.ones(3), angles=[40.0, 50.0, 60.0])
        assert np.array_equal(box.angles, [40.0, 50.0, 60.0])

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

    def test_angles_setter(self):
        box = mb.Box(mins=np.zeros(3), maxs=np.ones(3), angles=90*np.ones(3))
        box.angles = np.array([60.0, 120.0, 60.0])
        assert (box.angles == np.array([60.0, 120.0, 60.0])).all()

    def test_setters_with_lists(self):
        box = mb.Box(mins=np.zeros(3), maxs=2 * np.ones(3))
        box.mins = [1, 1, 1]
        assert (box.mins == np.ones(3)).all()
        assert (box.maxs - box.mins == np.ones(3)).all()
        assert (box.lengths == np.ones(3)).all()
        box.maxs = [3, 3, 3]
        assert (box.maxs == 3 * np.ones(3)).all()
        assert (box.maxs - box.mins == 2 * np.ones(3)).all()
        assert (box.lengths == 2 * np.ones(3)).all()
        box.lengths = [4, 4, 4]
        assert (box.lengths == 4 * np.ones(3)).all()
        assert (box.maxs - box.mins == 4 * np.ones(3)).all()
        box.angles = [90, 90, 120]
        assert (box.angles == np.array([90, 90, 120])).all()

    def test_single_dimension_setter(self):
        box = mb.Box(mins=np.zeros(3), maxs=4*np.ones(3))
        assert (box.lengths == 4*np.ones(3)).all()
        
        box.maxs[0] = 5
        box.mins[2] = 1
        assert np.allclose(box.mins, np.array([0, 0, 1], dtype=np.float))
        assert np.allclose(box.maxs, np.array([5, 4, 4], dtype=np.float))
        assert np.allclose(box.lengths, np.array([5, 4, 3], dtype=np.float))

        box.lengths[1] = 6
        assert np.allclose(box.mins, np.array([0, -1, 1], dtype=np.float))
        assert np.allclose(box.maxs, np.array([5, 5, 4], dtype=np.float))
        assert np.allclose(box.lengths, np.array([5, 6, 3], dtype=np.float))

        new_box = mb.Box(5)
        assert np.allclose(new_box.lengths, np.array([5, 5, 5], dtype=np.float))
        assert np.allclose(new_box.mins, np.array([0, 0, 0], dtype=np.float))
        assert np.allclose(new_box.maxs, np.array([5, 5, 5], dtype=np.float))

        new_box.lengths = 4
        assert np.allclose(new_box.lengths, np.array([4, 4, 4], dtype=np.float))
        assert np.allclose(new_box.mins, np.array([0.5, 0.5, 0.5], dtype=np.float))
        assert np.allclose(new_box.maxs, np.array([4.5, 4.5, 4.5], dtype=np.float))

    def test_sanity_checks(self):
        # Initialization step
        with pytest.raises(AssertionError):
            box = mb.Box(mins=[3,3,3], maxs=[1,1,1])
        with pytest.raises(AssertionError):
            box = mb.Box(lengths=-1)
        
        # Modifying step
        box = mb.Box(mins=[2,2,2], maxs=[4,4,4])
        with pytest.raises(AssertionError):
            box.mins[1] = box.maxs[1] + 1
        with pytest.raises(AssertionError): 
            box.maxs[1] = box.mins[1] - 1
        with pytest.raises(AssertionError):
            box.lengths = -1

    def test_compound_without_box(self, ethane):
        # Set coordinates to trigger the case where `box.mins`
        # coordinates can be less than [0, 0, 0]
        ethane.xyz = np.array([[0.31079999, 0.0653, -0.85259998],
                   [0.459, 0.0674, -0.8132],
                   [0.281, -0.0349, -0.8761],
                   [0.251, 0.1015, -0.7710],
                   [0.295, 0.1278, -0.9380],
                   [0.474, 0.0049, -0.7277],
                   [0.518, 0.0312, -0.8946],
                   [0.488, 0.1676, -0.7896]])
        # Set periodicity
        ethane.periodicity = np.array([0.767, 0.703, 0.710])
        # Convert compound to pmd.structure to set Box info
        # Conversion should happen without error
        ethane.to_parmed()
