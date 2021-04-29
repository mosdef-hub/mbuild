import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.exceptions import MBuildError


class TestBox(BaseTest):
    @pytest.mark.parametrize(
        "vectors",
        [
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[0.5, np.sqrt(3) / 2, 0.0], [0.5, -np.sqrt(3) / 2, 0], [0, 0, 1]],
        ],
    )
    def test_from_vectors(self, vectors):
        mb.Box.from_box_vectors(vectors)

    @pytest.mark.parametrize(
        "lengths, angles",
        [
            ([1, 1, 1], [90, 90, 90]),
            ([1, 5, 7], [90, 90, 90]),
            ([3, 4, 5], [90, 90, 120]),
        ],
    )
    def test_from_lengths_angles(self, lengths, angles):
        mb.Box(lengths=lengths, angles=angles)

    @pytest.mark.parametrize(
        "lh_matrix",
        [
            [[1, 2, 3], [3, 2, 1], [2, 1, 3]],
            [[1, 0, 0], [0, 1, 0], [0, -1, -1]],
            [[0.5, np.sqrt(3) / 2, 0.0], [0.5, 0.0, -np.sqrt(3) / 2], [0, 0, 1]],
        ],
    )
    def test_left_handed_matrix(self, lh_matrix):
        msg = (
            "Box vectors provided for a left-handed basis, these will "
            "be transformed into a right-handed basis automatically."
        )
        with pytest.warns(UserWarning, match=r"provided for a left\-handed basis"):
            # TODO add vector method properly
            mb.Box.from_box_vectors(vectors=lh_matrix)

    @pytest.mark.parametrize(
        "vecs",
        [
            [[1, 0, 0], [1, 1, 0], [1, 1, 0]],
            [[1, 0, 0], [5, 0, 0], [0, 0, 1]],
            [[1, 0, 0], [-5, 0, 0], [0, 0, 1]],
        ],
    )
    def test_colinear_vectors(self, vecs):
        with pytest.raises(
            mb.exceptions.MBuildError,
            match=r"co\-linear",
        ):
            mb.Box.from_box_vectors(vectors=vecs)

    @pytest.mark.parametrize(
        "vecs",
        [
            [[1.5, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[0.5, np.sqrt(3) / 2, 0.0], [0.5, -np.sqrt(3) / 2, 0], [0, 0, 8]],
        ],
    )
    def test_non_unit_vectors(self, vecs):
        lengths = [2, 2, 2]
        with pytest.raises(MBuildError, match=r"are not close to 1\.0"):
            mb.Box.from_uvec_lengths(uvec=vecs, lengths=lengths)
