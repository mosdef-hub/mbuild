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

    def test_target_coordinate_capture_radius_validation(self):
        # capture_radius needs a termination radius to ramp toward
        with pytest.raises(ValueError):
            TargetCoordinate(
                target_coordinate=(3, 3, 3), weight=1.0, capture_radius=2.0
            )
        # termination_radius is meaningless without capture_radius
        with pytest.raises(ValueError):
            TargetCoordinate(
                target_coordinate=(3, 3, 3), weight=1.0, termination_radius=0.3
            )
        # capture_radius must be outside termination_radius
        with pytest.raises(ValueError):
            TargetCoordinate(
                target_coordinate=(3, 3, 3),
                weight=1.0,
                capture_radius=0.3,
                termination_radius=0.5,
            )
        with pytest.raises(ValueError):
            TargetCoordinate(
                target_coordinate=(3, 3, 3),
                weight=1.0,
                capture_radius=-1.0,
                termination_radius=0.3,
            )
        with pytest.raises(ValueError):
            TargetCoordinate(
                target_coordinate=(3, 3, 3),
                weight=1.0,
                capture_radius="ten",
                termination_radius=0.3,
            )
        with pytest.raises(ValueError):
            TargetCoordinate(
                target_coordinate=(3, 3, 3),
                weight=0.5,
                capture_radius=2.0,
                termination_radius=0.3,
                adapt_sharpness=0.0,
            )

    def test_target_coordinate_terminator(self):
        terminator = WithinCoordinate(target_coordinate=(3, 3, 3), distance=0.3)
        bias = TargetCoordinate(
            target_coordinate=(3, 3, 3),
            weight=1.0,
            capture_radius=2.0,
            terminator=terminator,
        )
        assert bias.termination_radius == terminator.distance
        # terminator and termination_radius are mutually exclusive
        with pytest.raises(ValueError):
            TargetCoordinate(
                target_coordinate=(3, 3, 3),
                weight=1.0,
                capture_radius=2.0,
                termination_radius=0.3,
                terminator=terminator,
            )
        # terminator target must match the bias target
        with pytest.raises(ValueError):
            TargetCoordinate(
                target_coordinate=(1, 1, 1),
                weight=1.0,
                capture_radius=2.0,
                terminator=terminator,
            )

    def test_target_coordinate_auto_capture_radius(self):
        bias = TargetCoordinate(
            target_coordinate=(3, 3, 3),
            weight=1.0,
            capture_radius="auto",
            termination_radius=0.3,
        )
        assert bias._capture_radius is None
        hard_sphere_random_walk(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            bias=bias,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        assert bias._capture_radius == pytest.approx(2.5)
        # An auto radius that would fall inside the termination radius is caught
        bad_bias = TargetCoordinate(
            target_coordinate=(3, 3, 3),
            weight=1.0,
            capture_radius="auto",
            termination_radius=5.0,
        )
        with pytest.raises(ValueError):
            hard_sphere_random_walk(
                termination=Termination([NumSites(15), NumAttempts(1e4)]),
                bond_length=0.25,
                bias=bad_bias,
                initial_point=(0, 0, 0),
                radius=0.22,
                seed=14,
            )

    def test_target_coordinate_outside_capture_radius_matches_fixed_weight(self):
        """Outside capture_radius the walk is identical to a fixed-weight walk."""
        kwargs = dict(
            termination=Termination([NumSites(15), NumAttempts(1e4)]),
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        fixed = hard_sphere_random_walk(
            bias=TargetCoordinate(target_coordinate=(3, 3, 3), weight=0.8), **kwargs
        )
        # A capture_radius this walk never reaches leaves the ramp at r = 0
        never_captured = hard_sphere_random_walk(
            bias=TargetCoordinate(
                target_coordinate=(3, 3, 3),
                weight=0.8,
                capture_radius=0.301,
                termination_radius=0.3,
            ),
            **kwargs,
        )
        assert np.allclose(fixed.coordinates, never_captured.coordinates)

    def test_target_coordinate_capture_radius_ramp(self):
        """The effective weight ramps from weight at capture_radius up to 1."""
        bias = TargetCoordinate(
            target_coordinate=(0, 0, 0),
            weight=0.2,
            capture_radius=2.0,
            termination_radius=0.5,
        )
        assert bias.last_weight is None and bias.last_r is None

        def weight_at(distance):
            bias._adaptive_weight(np.array([distance, 0.0, 0.0]))
            return bias.last_r, bias.last_weight

        # Beyond capture_radius the ramp holds at weight
        assert weight_at(5.0) == (0.0, pytest.approx(0.2))
        assert weight_at(2.0) == (0.0, pytest.approx(0.2))
        # At and inside termination_radius the ramp is pinned at 1
        assert weight_at(0.5) == (1.0, pytest.approx(1.0))
        assert weight_at(0.0) == (1.0, pytest.approx(1.0))
        # Linear in between by default
        assert weight_at(1.25) == (pytest.approx(0.5), pytest.approx(0.6))
        # Monotonically increasing as the walk closes on the target
        weights = [weight_at(d)[1] for d in np.linspace(2.0, 0.5, 20)]
        assert np.all(np.diff(weights) >= 0)

    def test_target_coordinate_adapt_sharpness(self):
        """adapt_sharpness > 1 ramps later, < 1 ramps earlier."""

        def weight_at(sharpness, distance):
            bias = TargetCoordinate(
                target_coordinate=(0, 0, 0),
                weight=0.2,
                capture_radius=2.0,
                termination_radius=0.5,
                adapt_sharpness=sharpness,
            )
            bias._adaptive_weight(np.array([distance, 0.0, 0.0]))
            return bias.last_weight

        assert weight_at(3.0, 1.25) < weight_at(1.0, 1.25) < weight_at(0.3, 1.25)
        # Endpoints are unaffected by sharpness
        for sharpness in (0.3, 1.0, 3.0):
            assert weight_at(sharpness, 2.0) == pytest.approx(0.2)
            assert weight_at(sharpness, 0.5) == pytest.approx(1.0)

    def test_target_coordinate_capture_radius_reaches_target(self):
        """A ramped walk terminates inside the terminator's distance."""
        terminator = WithinCoordinate(target_coordinate=(2, 2, 2), distance=0.3)
        termination = Termination([terminator, NumAttempts(5e4)])
        kwargs = dict(
            bond_length=0.25,
            initial_point=(0, 0, 0),
            radius=0.22,
            seed=14,
        )
        adaptive = TargetCoordinate(
            target_coordinate=(2, 2, 2),
            weight=0.2,
            capture_radius=1.5,
            terminator=terminator,
        )
        path = hard_sphere_random_walk(termination=termination, bias=adaptive, **kwargs)
        end_distance = np.linalg.norm(path.coordinates[-1] - np.array([2, 2, 2]))
        assert end_distance <= 0.3 + terminator.tolerance

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

    def test_target_coordinate_capture_radius_pbc(self):
        """The capture-radius ramp measures distance to the target in minimum image."""
        from mbuild.path.path_utils import target_sq_distances

        L = 10.0
        target = np.array([-4.9, 0.0, 0.0])
        bias = TargetCoordinate(
            target_coordinate=target,
            weight=0.5,
            capture_radius=2.0,
            termination_radius=0.2,
        )

        class _State:
            pbc = np.array([True, True, True])
            box_lengths = np.array([L, L, L], dtype=np.float32)

        # 0.2 nm from the target across the +x boundary, 9.8 nm from it directly
        current = np.array([4.9, 0.0, 0.0])
        assert bias._adaptive_weight(current) == pytest.approx(0.5)

        bias.state = _State()
        bias._target_sq_distances = target_sq_distances
        assert bias._adaptive_weight(current) == pytest.approx(1.0)
