import pytest
import numpy as np
import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestBox(BaseTest):

    @pytest.mark.parametrize('vectors',
                                [
                                    [[1, 0, 0],[0, 1, 0],[0, 0, 1]],
                                    [[0.5, np.sqrt(3)/2, 0.0],
                                     [0.5, -np.sqrt(3)/2, 0],
                                     [0, 0, 1]],
                                ])
    def test_pass_working_box_vectors(self, vectors):
        print(mb.Box(box_vectors=vectors,))
        assert False

    @pytest.mark.parametrize('lengths, angles',
                             [
                                 ([1, 1, 1],[90, 90, 90]),
                                 ([1,5,7], [90, 90, 90]),
                                 ([3,4,5], [90, 90, 120])
                             ])
    def test_from_lengths_angles(self, lengths, angles):
        mb.Box.from_lengths_angles(lengths=lengths, angles=angles)

    @pytest.mark.parametrize('lh_matrix',
                             [
                                 [[1, 2, 3], [3, 2, 1], [2, 1, 3]]
                             ])
    def test_left_handed_matrix(self, lh_matrix):
        mb.Box(box_vectors=lh_matrix)

