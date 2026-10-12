"""Contains biases that can be included in mbuild.path.HardSphereRandomWalk."""

import numpy as np

from mbuild.path.path_utils import target_density as _target_density_cpu
from mbuild.path.path_utils import (
    target_sq_distances as _target_sq_distances_cpu,
)


class Bias:
    """Base class for biases applied to candidate sites in a random walk.

    Sub classes measure a signal over candidate sites and rank them with `_score`.

    Parameters
    ----------
    weight : float, required
        Bias weight in (0, 1]. 1 ranks candidates by the signal alone; values
        approaching 0 rank them close to randomly.
    """

    def __init__(self, weight):
        if weight <= 0 or weight > 1:
            raise ValueError(
                "weight should be larger than 0 and smaller than or equal to 1."
            )
        self.weight = weight
        self.beta, self.noise_scale = self._score_params(weight)

    @staticmethod
    def _score_params(weight):
        """Map a bias weight onto (beta, noise_scale) used to score candidates."""
        # Large beta diminishes the effect of noise
        return weight / max(1e-6, (1.0 - weight)), 1.0 - weight

    def _score(self, signal, weight=None):
        """Score candidates by their signal plus noise scaled to the signal's spread.

        Scores use `self.weight` unless `weight` overrides it for this step.
        A flat signal is ordered by noise alone.
        """
        if weight is None:
            beta, noise_scale = self.beta, self.noise_scale
        else:
            beta, noise_scale = self._score_params(weight)
        spread = np.std(signal)
        if spread == 0:
            spread = 1.0
        noise = self.rng.normal(0.0, noise_scale * spread, size=signal.shape)
        return beta * signal + noise

    def _attach_path(self, path, state):
        """Create access Path and RandomWalkState used by hard_sphere_random_walk."""
        self.path = path
        self.state = state
        # Inherit rng from the path for use in Bias classes
        self.rng = self.state.rng
        # Decide which path_utils implementations to use based on the path's device.
        self._target_sq_distances = _target_sq_distances_cpu
        self._target_density = _target_density_cpu

    def _clean(self):
        self.path = None
        self.state = None
        self.rng = None
        self._target_sq_distances = None
        self._target_density = None

    def __call__(self, candidates, coordinates, names):
        """Implemented in sub classes of Bias."""
        raise NotImplementedError


class TargetCoordinate(Bias):
    """Bias next-moves so that ones moving closer to a target coordinate are more likely to be accepted.

    Applies a fixed `weight` at every step unless `capture_radius` is set, in
    which case the weight at a distance `d` from the target is::

        r = clip((capture_radius - d) / (capture_radius - termination_radius), 0, 1)
        weight_eff = weight + (1 - weight) * r ** adapt_sharpness

    Parameters
    ----------
    target_coordinate : array-like (3,), required
        The target coordinate in units of nm.
    weight : float, required
        Bias weight in (0, 1]. Applied at every step when `capture_radius` is
        None, otherwise while outside `capture_radius`.
    capture_radius : float or "auto", optional, default=None
        Distance (nm) at which the weight begins ramping up. "auto" resolves to
        10 * bond_length when the bias is attached to a walk. Must be larger
        than `termination_radius`.
    termination_radius : float, optional, default=None
        Distance (nm) at which the weight reaches 1. Required when
        `capture_radius` is given, unless `terminator` supplies it.
    terminator : mbuild.path.termination.WithinCoordinate, optional, default=None
        Takes `termination_radius` from `terminator.distance` and validates
        `terminator.target_coordinate` against `target_coordinate`. Mutually
        exclusive with `termination_radius`.
    adapt_sharpness : float, optional, default=1.0
        Ramp exponent. 1.0 ramps linearly, > 1.0 ramps later and more sharply,
        < 1.0 ramps earlier and more gradually.

    Attributes
    ----------
    last_r : float or None
        Ramp fraction in [0, 1] used on the most recent step.
    last_weight : float or None
        Weight used on the most recent step.
    """

    def __init__(
        self,
        target_coordinate,
        weight,
        capture_radius=None,
        termination_radius=None,
        terminator=None,
        adapt_sharpness=1.0,
    ):
        self.target_coordinate = np.asarray(target_coordinate)
        self.last_r = None
        self.last_weight = None
        super(TargetCoordinate, self).__init__(weight=weight)

        if terminator is not None:
            if termination_radius is not None:
                raise ValueError(
                    "Pass only one of `termination_radius` and `terminator`; "
                    "`terminator.distance` is used as the termination radius."
                )
            if not np.allclose(terminator.target_coordinate, self.target_coordinate):
                raise ValueError(
                    "`terminator.target_coordinate` "
                    f"{np.asarray(terminator.target_coordinate)} does not match "
                    f"`target_coordinate` {self.target_coordinate}."
                )
            termination_radius = terminator.distance

        self.capture_radius = capture_radius
        self.termination_radius = termination_radius
        self.adapt_sharpness = adapt_sharpness
        # Resolved in _attach_path; equals capture_radius unless it is "auto".
        self._capture_radius = (
            None if isinstance(capture_radius, str) else capture_radius
        )

        if capture_radius is None:
            if termination_radius is not None:
                raise ValueError(
                    "`termination_radius`/`terminator` are only used when "
                    "`capture_radius` is set."
                )
            return

        if isinstance(capture_radius, str):
            if capture_radius != "auto":
                raise ValueError(
                    f"Unsupported {capture_radius=}; pass a float or 'auto'."
                )
        elif capture_radius <= 0:
            raise ValueError("`capture_radius` should be larger than 0.")

        if termination_radius is None:
            raise ValueError(
                "`capture_radius` requires `termination_radius` or `terminator` "
                "to define the end of the weight ramp."
            )
        if termination_radius < 0:
            raise ValueError("`termination_radius` should not be negative.")
        if adapt_sharpness <= 0:
            raise ValueError("`adapt_sharpness` should be larger than 0.")
        if not isinstance(capture_radius, str):
            self._validate_radii(capture_radius)

    def _validate_radii(self, capture_radius):
        if capture_radius <= self.termination_radius:
            raise ValueError(
                f"{capture_radius=} should be larger than "
                f"termination_radius={self.termination_radius}."
            )

    def _attach_path(self, path, state):
        """Create access Path and RandomWalkState used by hard_sphere_random_walk."""
        super(TargetCoordinate, self)._attach_path(path, state)
        if self.capture_radius == "auto":
            self._capture_radius = 10.0 * state.bond_length
            self._validate_radii(self._capture_radius)

    def _adaptive_weight(self, current_coordinate):
        """Return the weight for this step, ramped by distance to the target."""
        if getattr(self, "state", None) is None:
            distance = np.linalg.norm(self.target_coordinate - current_coordinate)
        else:
            sq_distance = self._target_sq_distances(
                self.target_coordinate,
                np.asarray(current_coordinate).reshape(1, 3),
                self.state.pbc,
                self.state.box_lengths,
            )
            distance = np.sqrt(sq_distance[0])
        span = max(self._capture_radius - self.termination_radius, 1e-12)
        r = float(np.clip((self._capture_radius - distance) / span, 0.0, 1.0))
        weight = self.weight + (1.0 - self.weight) * (r**self.adapt_sharpness)
        self.last_r = r
        self.last_weight = weight
        return weight

    def __call__(self, candidates, coordinates, names):
        """Sorts a set of candidate coordinates according to the bias.

        Parameters
        ----------
        candidates : np.ndarray (N,3), required
            Array of coordinates candidate sites for the next site in a random walk.
        coordinates : np.ndarray (N,3), required
            Array of coordinates for all previously accepted sites in a random walk.
        names : np.ndarray (N), required
            Array of site names for all previously accepted sites in a random walk.

        Returns
        -------
        candidates : np.ndarray (N,3)
            Returns the original candidate array, sorted according to the bias.
        """
        if self.capture_radius is None:
            weight = None
        else:
            weight = self._adaptive_weight(coordinates[-1])
        sq_distances = self._target_sq_distances(
            self.target_coordinate,
            candidates,
            self.state.pbc,
            self.state.box_lengths,
        )
        scores = self._score(sq_distances, weight=weight)
        # Target coordinate should favor short distances, sort in ascending order (np default)
        sort_idx = np.argsort(scores)
        return candidates[sort_idx]


class AvoidCoordinate(Bias):
    """Bias next-moves so that ones moving further from a specific coordinate are more likely to be accepted."""

    def __init__(self, avoid_coordinate, weight):
        self.avoid_coordinate = np.asarray(avoid_coordinate)
        super(AvoidCoordinate, self).__init__(weight=weight)

    def __call__(self, candidates, coordinates, names):
        """Sorts a set of candidate coordinates according to the bias.

        Parameters
        ----------
        candidates : np.ndarray (N,3), required
            Array of coordinates candidate sites for the next site in a random walk.
        coordinates : np.ndarray (N,3), required
            Array of coordinates for all previously accepted sites in a random walk.
        names : np.ndarray (N), required
            Array of site names for all previously accepted sites in a random walk.

        Returns
        -------
        candidates : np.ndarray (N,3)
            Returns the original candidate array, sorted according to the bias.
        """
        sq_distances = self._target_sq_distances(
            self.avoid_coordinate,
            candidates,
            self.state.pbc,
            self.state.box_lengths,
        )
        scores = self._score(sq_distances)
        # Avoid cooardinate should favor larger distances, sort in descending order
        sort_idx = np.argsort(scores)[::-1]
        return candidates[sort_idx]


class TargetType(Bias):
    """Bias next-moves so that ones increasing local density relative to a site type are more likely to be accepted."""

    def __init__(self, target_type, weight, r_cut):
        self.target_type = target_type
        self.r_cut = float(r_cut)
        super().__init__(weight=weight)

    def __call__(self, candidates, coordinates, names):
        """Sorts a set of candidate coordinates according to the bias.

        Parameters
        ----------
        candidates : np.ndarray (N,3), required
            Array of coordinates candidate sites for the next site in a random walk.
        coordinates : np.ndarray (N,3), required
            Array of coordinates for all previously accepted sites in a random walk.
        names : np.ndarray (N), required
            Array of site names for all previously accepted sites in a random walk.

        Returns
        -------
        candidates : np.ndarray (N,3)
            Returns the original candidate array, sorted according to the bias.
        """
        # Only get coordinates of sites where the site name == target type
        target_coords = coordinates[names == self.target_type]
        densities = self._target_density(
            candidates=candidates, target_coords=target_coords, r_cut=self.r_cut
        )
        scores = self._score(densities)
        # Target type should favor larger densities, sort in descending order
        sort_idx = np.argsort(scores)[::-1]
        return candidates[sort_idx]


class AvoidType(Bias):
    """Bias next-moves so that ones decreasing local density relative to a site type are more likely to be accepted."""

    def __init__(self, avoid_type, weight, r_cut):
        self.avoid_type = avoid_type
        self.r_cut = r_cut
        super().__init__(weight=weight)

    def __call__(self, candidates, coordinates, names):
        """Sorts a set of candidate coordinates according to the bias.

        Parameters
        ----------
        candidates : np.ndarray (N,3), required
            Array of coordinates candidate sites for the next site in a random walk.
        coordinates : np.ndarray (N,3), required
            Array of coordinates for all previously accepted sites in a random walk.
        names : np.ndarray (N), required
            Array of site names for all previously accepted sites in a random walk.

        Returns
        -------
        candidates : np.ndarray (N,3)
            Returns the original candidate array, sorted according to the bias.
        """
        target_coords = coordinates[names == self.avoid_type]
        densities = self._target_density(
            candidates=candidates, target_coords=target_coords, r_cut=self.r_cut
        )
        scores = self._score(densities)
        # Avoid type should favor smaller densities, sort in ascending order (np default)
        sort_idx = np.argsort(scores)
        return candidates[sort_idx]


class TargetDirection(Bias):
    """Bias next-moves so that ones moving along a target direction are more likely to be accepted."""

    def __init__(self, direction, weight):
        direction = np.asarray(direction, dtype=np.float32)
        norm = np.linalg.norm(direction)
        if norm < 1e-8:
            raise ValueError("Direction vector must be non-zero.")
        self.direction = direction / norm
        super().__init__(weight=weight)

    def __call__(self, candidates, coordinates, names):
        """Sorts a set of candidate coordinates according to the bias.

        Parameters
        ----------
        candidates : np.ndarray (N,3), required
            Array of coordinates candidate sites for the next site in a random walk.
        coordinates : np.ndarray (N,3), required
            Array of coordinates for all previously accepted sites in a random walk.
        names : np.ndarray (N), required
            Array of site names for all previously accepted sites in a random walk.

        Returns
        -------
        candidates : np.ndarray (N,3)
            Returns the original candidate array, sorted according to the bias.
        """
        last_step_pos = coordinates[-1]
        next_step_vectors = candidates - last_step_pos
        norms = np.linalg.norm(next_step_vectors, axis=1, keepdims=True)
        # Avoid divisions by zero
        norms = np.clip(norms, 1e-12, None)
        next_step_unit_vectors = next_step_vectors / norms
        # Alignment score: dot product with target direction
        alignment = np.dot(next_step_unit_vectors, self.direction)
        scores = self._score(alignment)
        # Larger dot product = better alignment with target, sort descending
        sort_idx = np.argsort(scores)[::-1]
        return candidates[sort_idx]


class AvoidDirection(Bias):
    """Bias next-moves so that ones not moving along a certain direction are more likely to be accepted."""

    def __init__(self, direction, weight):
        self.direction = direction
        super().__init__(weight=weight)

    def __call__(self, candidates, coordinates, names):
        """Sorts a set of candidate coordinates according to the bias.

        Parameters
        ----------
        candidates : np.ndarray (N,3), required
            Array of coordinates candidate sites for the next site in a random walk.
        coordinates : np.ndarray (N,3), required
            Array of coordinates for all previously accepted sites in a random walk.
        names : np.ndarray (N), required
            Array of site names for all previously accepted sites in a random walk.

        Returns
        -------
        candidates : np.ndarray (N,3)
            Returns the original candidate array, sorted according to the bias.
        """
        last_step_pos = coordinates[-1]
        next_step_vectors = candidates - last_step_pos
        norms = np.linalg.norm(next_step_vectors, axis=1, keepdims=True)
        # Avoid divisions by zero
        norms = np.clip(norms, 1e-12, None)
        next_step_unit_vectors = next_step_vectors / norms
        # Alignment score: dot product with target direction
        alignment = np.dot(next_step_unit_vectors, self.direction)
        scores = self._score(alignment)
        # Larger dot product = better alignment with target, sort ascending (np default)
        sort_idx = np.argsort(scores)
        return candidates[sort_idx]
