import numpy as np

from mbuild.path import HardSphereRandomWalk
from mbuild.path.bias import (
    AvoidCoordinate,
    AvoidDirection,
    AvoidType,
    TargetCoordinate,
    TargetDirection,
    TargetType,
)
from mbuild.tests.base_test import BaseTest


def radius_of_gyration(coordinates):
    """Calculate the radius of gyration for a set of coordinates using the geometric center."""
    coordinates = np.array(coordinates)
    geometric_center = np.mean(coordinates, axis=0)
    squared_distances = np.sum((coordinates - geometric_center) ** 2, axis=1)
    rg = np.sqrt(np.mean(squared_distances))
    return rg


class TestBias(BaseTest):
    def test_target_coordinate(self):
        bias = TargetCoordinate(target_coordinate=(3, 3, 3), weight=1.0)
        rw_path = HardSphereRandomWalk(
            N=15,
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        dist_to_target = np.linalg.norm(rw_path.coordinates[-1] - np.array([3, 3, 3]))
        rw_path_biased = HardSphereRandomWalk(
            N=15,
            bond_length=0.25,
            bias=bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        dist_to_target_biased = np.linalg.norm(
            rw_path_biased.coordinates[-1] - np.array([3, 3, 3])
        )
        assert dist_to_target_biased < dist_to_target

    def test_avoid_coordinate(self):
        bias = AvoidCoordinate(avoid_coordinate=(3, 3, 3), weight=1.0)
        rw_path = HardSphereRandomWalk(
            N=15,
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        dist_to_target = np.linalg.norm(rw_path.coordinates[-1] - np.array([3, 3, 3]))
        rw_path_biased = HardSphereRandomWalk(
            N=15,
            bond_length=0.25,
            bias=bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
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
            N=50,
            bias=target_bias,
            bond_length=0.25,
            initial_point=(0, 0, 0),
            bead_name="A",
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        rw_path_avoid = HardSphereRandomWalk(
            N=50,
            bond_length=0.25,
            bias=avoid_bias,
            initial_point=(0, 0, 0),
            bead_name="A",
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        assert radius_of_gyration(rw_path_target.coordinates) < radius_of_gyration(
            rw_path_avoid.coordinates
        )

    def test_target_direction(self):
        target_bias = TargetDirection(direction=(1, 0, 0), weight=0.7)
        avoid_bias = AvoidDirection(direction=(1, 0, 0), weight=0.7)
        rw_path = HardSphereRandomWalk(
            N=15,
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        head_tail_vec = rw_path.coordinates[-1] - rw_path.coordinates[0]

        rw_path_target = HardSphereRandomWalk(
            N=15,
            bond_length=0.25,
            bias=target_bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        head_tail_vec_target = (
            rw_path_target.coordinates[-1] - rw_path_target.coordinates[0]
        )

        rw_path_avoid = HardSphereRandomWalk(
            N=15,
            bond_length=0.25,
            bias=avoid_bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
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
