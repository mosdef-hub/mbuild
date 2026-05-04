import time

import numpy as np
import pytest

from mbuild.path.bias import (
    TargetCoordinate,
)
from mbuild.path.build import hard_sphere_random_walk
from mbuild.path.constraints import (
    CuboidConstraint,
)
from mbuild.path.termination import (
    ContourLength,
    EndToEndDistance,
    NumAttempts,
    NumSites,
    RadiusOfGyration,
    Termination,
    Terminator,
    WallTime,
    WithinCoordinate,
)
from mbuild.tests.base_test import BaseTest, radius_of_gyration


class TestTermination(BaseTest):
    def test_no_termination(self):
        with pytest.raises(ValueError):
            hard_sphere_random_walk(
                termination=None,
                bond_length=0.25,
                radius=0.22,
                seed=14,
                chunk_size=50,
            )

    def test_terminator_base_class(self):
        with pytest.raises(NotImplementedError):
            termination = Termination(Terminator(is_target=True))
            termination.is_met(coordinates=None, names=None)

    def test_single_terminator(self):
        num_sites = NumSites(5000)
        termination = Termination(num_sites)
        assert isinstance(termination.terminators, list)
        assert termination.terminators[0] is num_sites

    def test_contour_length_termination(self):
        rw = hard_sphere_random_walk(
            radius=0.20,
            bond_length=0.22,
            termination=Termination([ContourLength(22), NumAttempts(1e4)]),
            initial_point=(0, 0, 0),
        )
        assert len(rw.coordinates) == 100

    def test_wall_time_termination(self):
        cube = CuboidConstraint(Lx=1, Ly=1, Lz=1)
        start = time.time()
        rw = hard_sphere_random_walk(
            radius=0.20,
            bond_length=0.22,
            volume_constraint=cube,
            termination=Termination([NumSites(5000), WallTime(5)]),
            initial_point=(0, 0, 0),
            seed=14,
        )
        # Allow some buffer
        assert time.time() - start <= 6
        assert len(rw.coordinates) == 0

    def test_num_attempts_termination(self):
        rw_path = hard_sphere_random_walk(
            termination=Termination([NumSites(500), NumAttempts(100)]),
            bond_length=0.25,
            radius=0.22,
            seed=14,
        )
        assert len(rw_path.coordinates) == 0

    def test_rg_termination(self):
        num_sites = NumSites(20)
        rg = RadiusOfGyration(2)
        max_attempts = NumAttempts(5e4)
        termination = Termination([num_sites, rg, max_attempts])
        rw = hard_sphere_random_walk(
            termination=termination,
            initial_point=(0, 0, 0),
            bond_length=0.24,
            radius=0.22,
            seed=10,
        )
        assert len(rw.coordinates) >= 20
        assert np.allclose(radius_of_gyration(rw.coordinates), 2, atol=1e-1)

    def test_re_termination(self):
        num_sites = NumSites(20)
        re = EndToEndDistance(2, tolerance=0.1)
        max_attempts = NumAttempts(2e4)
        termination = Termination([num_sites, re, max_attempts])
        rw = hard_sphere_random_walk(
            termination=termination,
            bond_length=0.25,
            radius=0.22,
            seed=14,
        )
        dist = np.linalg.norm(rw.coordinates[-1] - rw.coordinates[0])
        assert np.allclose(dist, 2, atol=0.1)
        assert len(rw.coordinates) >= 20

    def test_within_coord_termination(self):
        bias = TargetCoordinate(target_coordinate=(2, 2, 2), weight=0.8)
        termination = Termination(
            [
                WithinCoordinate(target_coordinate=(2, 2, 2), distance=0.22),
                NumAttempts(1e3),
            ]
        )
        rw = hard_sphere_random_walk(
            termination=termination,
            bias=bias,
            initial_point=(0, 0, 0),
            bond_length=0.25,
            radius=0.22,
            trial_batch_size=100,
            seed=14,
        )
        dist = np.linalg.norm(rw.coordinates[-1] - np.array([2, 2, 2]))
        assert dist <= 0.22

        termination = Termination(
            [
                WithinCoordinate(
                    target_coordinate=(3, 3, 3), distance=0.0, tolerance=1e-1
                ),
                NumAttempts(100),
            ]
        )
        bias = TargetCoordinate(target_coordinate=(3, 3, 3), weight=1.0)
        rw = hard_sphere_random_walk(
            termination=termination,
            bias=bias,
            initial_point=(0, 0, 0),
            bond_length=0.25,
            radius=0.22,
            trial_batch_size=200,
            seed=14,
        )
        dist = np.linalg.norm(rw.coordinates[-1] - np.array([3, 3, 3]))
        assert np.allclose(dist, 0, atol=1e-1)
