import time

import numpy as np
import pytest

from mbuild.path import (
    HardSphereRandomWalk,
)
from mbuild.path.bias import TargetCoordinate
from mbuild.path.termination import (
    EndToEndDistance,
    NumAttempts,
    NumSites,
    RadiusOfGyration,
    Termination,
    WallTime,
    WithinCoordinate,
)
from mbuild.tests.base_test import BaseTest, radius_of_gyration
from mbuild.utils.volumes import (
    CuboidConstraint,
)


class TestTermination(BaseTest):
    def test_no_termination(self):
        with pytest.raises(RuntimeError):
            HardSphereRandomWalk(
                termination=None,
                bond_length=0.25,
                radius=0.22,
                min_angle=np.pi / 4,
                max_angle=np.pi,
                max_attempts=1e4,
                seed=14,
            )

    def test_wall_time_termination(self):
        cube = CuboidConstraint(Lx=1, Ly=1, Lz=1)
        start = time.time()
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(5000), WallTime(5)]),
            bond_length=0.25,
            radius=0.22,
            volume_constraint=cube,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        # Allow some buffer
        assert time.time() - start <= 6
        assert len(rw_path.coordinates) < 5000

    def test_num_attempts_termination(self):
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(500), NumAttempts(100)]),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        # First 2 steps in unconstrianed RWs are always accepted
        # Counting doesn't start until afterwards
        assert len(rw_path.coordinates) == 102

    def test_rg_termination(self):
        num_sites = NumSites(20)
        rg = RadiusOfGyration(3)
        max_attempts = NumAttempts(1e4)
        termination = Termination([num_sites, rg, max_attempts])
        rw_path = HardSphereRandomWalk(
            termination=termination,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        assert np.allclose(radius_of_gyration(rw_path.coordinates), 3, atol=1e-1)

    def test_re_termination(self):
        num_sites = NumSites(20)
        re = EndToEndDistance(2)
        max_attempts = NumAttempts(1e4)
        termination = Termination([num_sites, re, max_attempts])
        rw_path = HardSphereRandomWalk(
            termination=termination,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        dist = np.linalg.norm(rw_path.coordinates[-1] - rw_path.coordinates[0])
        assert np.allclose(dist, 2, atol=1e-1)

    def test_within_coord_termination(self):
        bias = TargetCoordinate(target_coordinate=(2, 2, 2), weight=0.8)
        termination = Termination(
            [
                WithinCoordinate(target_coordinate=(2, 2, 2), distance=0.22),
                NumAttempts(1e3),
            ]
        )
        rw_path = HardSphereRandomWalk(
            termination=termination,
            bias=bias,
            initial_point=(0, 0, 0),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            trial_batch_size=100,
            seed=14,
        )
        dist = np.linalg.norm(rw_path.coordinates[-1] - np.array([2, 2, 2]))
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
        rw_path = HardSphereRandomWalk(
            termination=termination,
            bias=bias,
            initial_point=(0, 0, 0),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=5e4,
            trial_batch_size=200,
            seed=14,
        )
        dist = np.linalg.norm(rw_path.coordinates[-1] - np.array([3, 3, 3]))
        assert np.allclose(dist, 0, atol=1e-1)
