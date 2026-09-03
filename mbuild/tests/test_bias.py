import numpy as np
import pytest

from mbuild.path.bias import (
    AvoidCoordinate,
    AvoidDirection,
    AvoidType,
    TargetCoordinate,
    TargetDirection,
    TargetType,
)
from mbuild.path.build import hard_sphere_random_walk
from mbuild.path.constraints import CuboidConstraint
from mbuild.path.termination import (
    NumAttempts,
    NumSites,
    Termination,
    WithinCoordinate,
)
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
        rw_path = hard_sphere_random_walk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        dist_to_target = np.linalg.norm(rw_path.coordinates[-1] - np.array([3, 3, 3]))
        rw_path_biased = hard_sphere_random_walk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        dist_to_target_biased = np.linalg.norm(
            rw_path_biased.coordinates[-1] - np.array([3, 3, 3])
        )
        assert dist_to_target_biased < dist_to_target

    def test_avoid_coordinate(self):
        bias = AvoidCoordinate(avoid_coordinate=(3, 3, 3), weight=1.0)
        rw_path = hard_sphere_random_walk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        dist_to_target = np.linalg.norm(rw_path.coordinates[-1] - np.array([3, 3, 3]))
        rw_path_biased = hard_sphere_random_walk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        dist_to_target_biased = np.linalg.norm(
            rw_path_biased.coordinates[-1] - np.array([3, 3, 3])
        )
        assert dist_to_target_biased > dist_to_target

    def test_target_and_avoid_type(self):
        target_bias = TargetType(target_type="A", weight=0.6, r_cut=2)
        avoid_bias = AvoidType(avoid_type="A", weight=0.6, r_cut=2)

        rw_path_target = hard_sphere_random_walk(
            termination=Termination([NumSites(50), NumAttempts(1e4)]),
            bias=target_bias,
            bond_length=0.25,
            initial_point=(0, 0, 0),
            bead_name="A",
            radius=0.22,
            seed=14,
        )
        rw_path_avoid = hard_sphere_random_walk(
            termination=Termination([NumSites(50), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=avoid_bias,
            initial_point=(0, 0, 0),
            bead_name="A",
            radius=0.22,
            seed=14,
        )
        assert radius_of_gyration(rw_path_target.coordinates) < radius_of_gyration(
            rw_path_avoid.coordinates
        )

    def test_target_direction(self):
        target_bias = TargetDirection(direction=(1, 0, 0), weight=0.7)
        avoid_bias = AvoidDirection(direction=(1, 0, 0), weight=0.7)
        rw_path = hard_sphere_random_walk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        head_tail_vec = rw_path.coordinates[-1] - rw_path.coordinates[0]

        rw_path_target = hard_sphere_random_walk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=target_bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        head_tail_vec_target = (
            rw_path_target.coordinates[-1] - rw_path_target.coordinates[0]
        )

        rw_path_avoid = hard_sphere_random_walk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=avoid_bias,
            initial_point=(0, 0, 0),
            radius=0.22,
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

    def test_score_affine_invariance(self):
        """Scores rank candidates identically under any positive affine rescale."""
        bias = TargetCoordinate(target_coordinate=(3, 3, 3), weight=0.6)
        signal = np.array([0.1, 4.2, 1.7, 3.3, 0.9, 2.5])
        for scale, offset in [(1.0, 0.0), (1e3, 0.0), (1e-3, 0.0), (7.0, 250.0)]:
            bias.rng = np.random.default_rng(5)
            baseline = np.argsort(bias._score(signal))
            bias.rng = np.random.default_rng(5)
            rescaled = np.argsort(bias._score(signal * scale + offset))
            assert np.array_equal(baseline, rescaled)

    def test_score_flat_signal_is_random(self):
        """A signal carrying no information leaves ordering to noise."""
        bias = TargetCoordinate(target_coordinate=(3, 3, 3), weight=0.6)
        bias.rng = np.random.default_rng(5)
        orders = {tuple(np.argsort(bias._score(np.full(6, 2.0)))) for _ in range(50)}
        assert len(orders) > 1

    def test_score_weight_controls_signal_to_noise(self):
        """Higher weight follows the signal more closely, at any signal scale."""
        signal = np.arange(20, dtype=float)
        best = np.argmin(signal)
        for scale in (1e-3, 1.0, 1e3):
            picked = {}
            for weight in (0.2, 0.9):
                bias = TargetCoordinate(target_coordinate=(0, 0, 0), weight=weight)
                bias.rng = np.random.default_rng(5)
                picked[weight] = np.mean(
                    [
                        np.argsort(bias._score(signal * scale))[0] == best
                        for _ in range(400)
                    ]
                )
            assert picked[0.9] > picked[0.2]
            # A given weight behaves the same at every signal scale
            assert picked[0.9] > 0.85
            assert 0.1 < picked[0.2] < 0.6

    def test_target_coordinate_pbc(self):
        """Under PBC the bias steers across the nearest boundary, not through the box."""
        L = 10.0
        box = CuboidConstraint(
            Lx=L, Ly=L, Lz=L, center=(0, 0, 0), pbc=(True, True, True)
        )
        start = np.array([4.0, 0.0, 0.0])
        target = np.array([-4.0, 0.0, 0.0])
        # target is 8 nm away directly, but only 2 nm across the +x boundary
        path = hard_sphere_random_walk(
            termination=Termination([NumSites(6), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.2,
            bias=TargetCoordinate(target_coordinate=target, weight=0.95),
            initial_point=start,
            volume_constraint=box,
            seed=3,
        )
        # Steps head toward the near boundary (+x), away from the raw target
        first_steps = path.coordinates[1:5, 0] - path.coordinates[0:4, 0]
        assert np.all(first_steps > 0)

    def test_within_coordinate_pbc(self):
        """WithinCoordinate measures distance to the target in minimum image."""
        L = 10.0
        target = np.array([-4.9, 0.0, 0.0])
        # 0.2 nm from the target across the +x boundary, 9.8 nm from it directly
        coordinates = np.array([[4.9, 0.0, 0.0]])
        names = np.array(["_A"])

        class _State:
            pbc = np.array([True, True, True])
            box_lengths = np.array([L, L, L], dtype=np.float32)

        unattached = WithinCoordinate(target_coordinate=target, distance=0.3)
        assert not unattached.is_met(coordinates, names)

        periodic = WithinCoordinate(target_coordinate=target, distance=0.3)
        periodic.state = _State()
        assert periodic.is_met(coordinates, names)
