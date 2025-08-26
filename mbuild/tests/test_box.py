import logging

import numpy as np
import pytest

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.tests.base_test import BaseTest


class TestBox(BaseTest):
    @pytest.mark.parametrize(
        "vectors",
        [
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[0.5, np.sqrt(3) / 2, 0.0], [0.5, -np.sqrt(3) / 2, 0], [0, 0, 1]],
        ],
    )
    def test_from_vectors(self, vectors):
        mb.Box.from_vectors(vectors)

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
        mb.Box.from_lengths_angles(lengths=lengths, angles=angles)

    @pytest.mark.parametrize(
        "lh_matrix",
        [
            [[1, 2, 3], [3, 2, 1], [2, 1, 3]],
            [[1, 0, 0], [0, 1, 0], [0, -1, -1]],
            [
                [0.5, np.sqrt(3) / 2, 0.0],
                [0.5, 0.0, -np.sqrt(3) / 2],
                [0, 0, 1],
            ],
        ],
    )
    def test_left_handed_matrix(self, lh_matrix, caplog):
        with caplog.at_level(logging.WARNING, logger="mbuild"):
            mb.Box.from_vectors(vectors=lh_matrix)
        assert "provided for a left-handed basis" in caplog.text

    @pytest.mark.parametrize(
        "vecs",
        [
            [[1, 0, 0], [1, 1, 0], [1, 1, 0]],
            [[1, 0, 0], [5, 0, 0], [0, 0, 1]],
            [[1, 0, 0], [-5, 0, 0], [0, 0, 1]],
        ],
    )
    def test_colinear_vectors(self, vecs):
        with pytest.raises(mb.exceptions.MBuildError, match=r"co\-linear"):
            mb.Box.from_vectors(vectors=vecs)

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

    @pytest.mark.parametrize(
        "uvec, lengths, expected_angles",
        [
            ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [3, 3, 3], [90, 90, 90]),
            (
                [
                    [np.sqrt(3) / 2, -0.5, 0],
                    [np.sqrt(3) / 2, 0.5, 0],
                    [0, 0, 1],
                ],
                [3, 3, 3],
                [90, 90, 60],
            ),
        ],
    )
    def test_uvec_lengths(self, uvec, lengths, expected_angles):
        box = mb.Box.from_uvec_lengths(uvec, lengths)
        assert np.all(np.isclose(box.lengths, lengths))
        assert np.all(np.isclose(box.angles, expected_angles))

    @pytest.mark.parametrize(
        "lengths, tilt_factors, angles",
        [
            ([1, 1, 1], [0.0, 0.0, 0.0], [90, 90, 90]),
            ([1.0, 1.0, 1.0], [-0.57735, 0.0, 0.0], [90.0, 90.0, 120.0]),
            ([3.0, 3.0, 1.0], [-0.57735, 0.0, 0.0], [90.0, 90.0, 120.0]),
        ],
    )
    def test_lengths_tilt_factors(self, lengths, tilt_factors, angles):
        box = mb.Box.from_lengths_tilt_factors(
            lengths,
            tilt_factors,
        )
        assert np.all(np.isclose(box.lengths, lengths))
        assert np.all(np.isclose(box.tilt_factors, tilt_factors))
        assert np.all(np.isclose(box.angles, angles))

    @pytest.mark.parametrize(
        "lo, hi, tilt_factors, angles",
        [
            ([0, 0, 0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [90, 90, 90]),
            (
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [-0.57735, 0.0, 0.0],
                [90.0, 90.0, 120.0],
            ),
            (
                [1.0, 2.0, -1.0],
                [2.0, 4.0, 3.0],
                [-0.57735, 0.0, 0.0],
                [90.0, 90.0, 120.0],
            ),
        ],
    )
    def test_lo_hi_tilt_factors(self, lo, hi, tilt_factors, angles):
        (xlo, ylo, zlo) = lo
        (xhi, yhi, zhi) = hi
        (xy, xz, yz) = tilt_factors

        lengths = [xhi - (xlo + xy), yhi - ylo, zhi - zlo]
        box = mb.Box.from_lo_hi_tilt_factors(lo=lo, hi=hi, tilt_factors=tilt_factors)
        assert np.all(np.isclose(box.lengths, lengths))
        assert np.all(np.isclose(box.tilt_factors, tilt_factors))
        assert np.all(np.isclose(box.angles, angles))

    @pytest.mark.parametrize(
        "a, b, c, alpha, beta, gamma",
        [(1, 1, 1, 90, 90, 90), (1, 1, 1, 90, 90, 120), (3, 6, 7, 97, 99, 120)],
    )
    def test_bravais_parameters(self, a, b, c, alpha, beta, gamma):
        box = mb.Box(lengths=[a, b, c], angles=[alpha, beta, gamma])
        assert np.all(
            np.isclose(list(box.bravais_parameters), [a, b, c, alpha, beta, gamma])
        )
