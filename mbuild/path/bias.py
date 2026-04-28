"""Contains biases that can be included in mbuild.path.HardSphereRandomWalk."""

import numpy as np

from mbuild.path.path_utils import target_density as _target_density_cpu
from mbuild.path.path_utils import (
    target_sq_distances as _target_sq_distances_cpu,
)

try:
    from numba import cuda

    _CUDA_AVAILABLE = cuda.is_available()
    if _CUDA_AVAILABLE:
        from mbuild.path.path_utils_gpu import (
            target_density as _target_density_gpu,
        )
        from mbuild.path.path_utils_gpu import (
            target_sq_distances as _target_sq_distances_gpu,
        )
    else:
        _target_density_gpu = None
        _target_sq_distances_gpu = None

except Exception:  # pragma: no cover - CUDA stack not importable or GPU utils failed
    _CUDA_AVAILABLE = False
    _target_density_gpu = None
    _target_sq_distances_gpu = None


class Bias:
    def __init__(self, weight):
        if weight <= 0 or weight > 1:
            raise ValueError(
                "weight should be larger than 0 and smaller than or equal to 1."
            )
        self.weight = weight
        # Large beta diminishes the effect of noise
        self.beta = self.weight / max(1e-6, (1.0 - self.weight))
        self.noise_scale = 1 - self.weight

    def _attach_path(self, path, state):
        """Create access Path and RandomWalkState used by hard_sphere_random_walk."""
        self.path = path
        self.state = state
        # Inherit rng from the path for use in Bias classes
        self.rng = self.state.rng
        # Decide which path_utils implementations to use based on the path's device.
        use_gpu = getattr(state, "run_on_gpu", False)
        if use_gpu and _CUDA_AVAILABLE and _target_sq_distances_gpu is not None:
            self._target_sq_distances = _target_sq_distances_gpu
            self._target_density = _target_density_gpu
        else:
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
    """Bias next-moves so that ones moving closer to a target coordinate are more likely to be accepted."""

    def __init__(self, target_coordinate, weight):
        self.target_coordinate = np.asarray(target_coordinate)
        super(TargetCoordinate, self).__init__(weight=weight)

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
        sq_distances = self._target_sq_distances(self.target_coordinate, candidates)
        noise = self.rng.normal(0, self.noise_scale, size=sq_distances.shape)
        scores = self.beta * sq_distances + noise
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
        sq_distances = self._target_sq_distances(self.avoid_coordinate, candidates)
        noise = self.rng.normal(0, self.noise_scale, size=sq_distances.shape)
        scores = self.beta * sq_distances + noise
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
        noise = self.rng.normal(0, self.noise_scale, size=densities.shape)
        scores = self.beta * densities + noise
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
        noise = self.rng.normal(0, self.noise_scale, size=densities.shape)
        scores = self.beta * densities + noise
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
        noise = self.rng.normal(0.0, self.noise_scale, size=alignment.shape)
        scores = self.beta * alignment + noise
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
        noise = self.rng.normal(0.0, self.noise_scale, size=alignment.shape)
        scores = self.beta * alignment + noise
        # Larger dot product = better alignment with target, sort ascending (np default)
        sort_idx = np.argsort(scores)
        return candidates[sort_idx]
