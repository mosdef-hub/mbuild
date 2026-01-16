import numpy as np
import pytest

from mbuild.path import HardSphereRandomWalk
from mbuild.path.bias import (
    AvoidCoordinate,
    AvoidDirection,
    AvoidType,
    TargetCoordinate,
    TargetDirection,
    TargetType,
)
from mbuild.path.termination import NumAttempts, NumSites, Termination
from mbuild.tests.base_test import BaseTest, radius_of_gyration


class TestBias(BaseTest):
    def test_bad_weight(self):
        with pytest.raises(ValueError):
            AvoidType(avoid_type="A", weight=0.0, r_cut=1)

        with pytest.raises(ValueError):
            AvoidType(avoid_type="A", weight=2.0, r_cut=1)

        with pytest.raises(ValueError):
            AvoidType(avoid_type="A", weight=-0.5, r_cut=1)

    def test_target_coordinate(self):
        bias = TargetCoordinate(target_coordinate=(3, 3, 3), weight=1.0)
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        dist_to_target = np.linalg.norm(rw_path.coordinates[-1] - np.array([3, 3, 3]))
        rw_path_biased = HardSphereRandomWalk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        dist_to_target_biased = np.linalg.norm(
            rw_path_biased.coordinates[-1] - np.array([3, 3, 3])
        )
        assert dist_to_target_biased < dist_to_target

    def test_avoid_coordinate(self):
        bias = AvoidCoordinate(avoid_coordinate=(3, 3, 3), weight=1.0)
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        dist_to_target = np.linalg.norm(rw_path.coordinates[-1] - np.array([3, 3, 3]))
        rw_path_biased = HardSphereRandomWalk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        dist_to_target_biased = np.linalg.norm(
            rw_path_biased.coordinates[-1] - np.array([3, 3, 3])
        )
        assert dist_to_target_biased > dist_to_target

    def test_target_and_avoid_type(self):
        target_bias = TargetType(target_type="A", weight=0.6, r_cut=2)
        avoid_bias = AvoidType(avoid_type="A", weight=0.6, r_cut=2)

        rw_path_target = HardSphereRandomWalk(
            termination=Termination([NumSites(50), NumAttempts(1e4)]),
            bias=target_bias,
            bond_length=0.25,
            initial_point=(0, 0, 0),
            bead_name="A",
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        rw_path_avoid = HardSphereRandomWalk(
            termination=Termination([NumSites(50), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=avoid_bias,
            initial_point=(0, 0, 0),
            bead_name="A",
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        assert radius_of_gyration(rw_path_target.coordinates) < radius_of_gyration(
            rw_path_avoid.coordinates
        )

    def test_target_direction(self):
        target_bias = TargetDirection(direction=(1, 0, 0), weight=0.7)
        avoid_bias = AvoidDirection(direction=(1, 0, 0), weight=0.7)
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        head_tail_vec = rw_path.coordinates[-1] - rw_path.coordinates[0]

        rw_path_target = HardSphereRandomWalk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=target_bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        head_tail_vec_target = (
            rw_path_target.coordinates[-1] - rw_path_target.coordinates[0]
        )

        rw_path_avoid = HardSphereRandomWalk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=avoid_bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            seed=14,
        )
        head_tail_vec_avoid = (
            rw_path_avoid.coordinates[-1] - rw_path_avoid.coordinates[0]
        )

        assert np.dot(head_tail_vec_target, np.array([1, 0, 0])) > np.dot(
            head_tail_vec, np.array([1, 0, 0])
        )
        assert np.dot(head_tail_vec_avoid, np.array([1, 0, 0])) < np.dot(
            head_tail_vec, np.array([1, 0, 0])
        )
        assert np.dot(head_tail_vec_avoid, np.array([1, 0, 0])) < np.dot(
            head_tail_vec_target, np.array([1, 0, 0])
        )
