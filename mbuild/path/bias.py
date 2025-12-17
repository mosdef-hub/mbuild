"""Contains biases that can be included in mbuild.path.HardSphereRandomWalk."""

import numpy as np

from mbuild.path.path_utils import target_density, target_sq_distances


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

    def _attach_path(self, path):
        self.path = path
        # Inherit rng from the path for use in Bias classes
        self.rng = self.path.rng

    def __call__(self, candidates):
        """Implemented in sub classes of Bias."""
        raise NotImplementedError


class TargetCoordinate(Bias):
    """Bias next-moves so that ones moving closer to a target coordinate are more likely to be accepted."""

    def __init__(self, target_coordinate, weight):
        self.target_coordinate = np.asarray(target_coordinate)
        super(TargetCoordinate, self).__init__(weight=weight)

    def __call__(self, candidates):
        sq_distances = target_sq_distances(self.target_coordinate, candidates)
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

    def __call__(self, candidates):
        sq_distances = target_sq_distances(self.avoid_coordinate, candidates)
        noise = self.rng.normal(0, self.noise_scale, size=sq_distances.shape)
        scores = self.beta * sq_distances + noise
        # Avoid cooardinate should favor larger distances, sort in descending order
        sort_idx = np.argsort(scores)[::-1]
        return candidates[sort_idx]


class TargetType(Bias):
    """Bias next-moves so that ones increasing local density relative to a site type are more likely to be accepted."""

    def __init__(self, target_type, weight, r_cut):
        self.target_type = target_type
        self.r_cut = r_cut
        super(TargetType, self).__init__(weight=weight)

    def __call__(self, candidates):
        types = np.array(
            [node[1]["name"] for node in self.path.bond_graph.nodes(data=True)]
        )
        # path.coordinates is an array with place holder values for future sites
        # Only check as far as types are defined (i.e., once a node is added from an accepted move.)
        target_coords = self.path.coordinates[: len(types)][
            types == self.target_type
        ].astype(np.float32)
        densities = target_density(
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
        super(AvoidType, self).__init__(weight=weight)

    def __call__(self, candidates):
        types = np.array(
            [node[1]["name"] for node in self.path.bond_graph.nodes(data=True)]
        )
        # path.coordinates is an array with place holder values for future sites
        # Only check as far as types are defined (i.e., once a node is added)
        target_coords = self.path.coordinates[: len(types)][
            types == self.avoid_type
        ].astype(np.float32)
        densities = target_density(
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
        super(TargetDirection, self).__init__(weight=weight)

    def __call__(self, candidates):
        last_step_pos = self.path.coordinates[self.path.count]
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
        super(AvoidDirection, self).__init__(weight=weight)

    def __call__(self, candidates):
        last_step_pos = self.path.coordinates[self.path.count]
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


class TargetEdge(Bias):
    """Bias next-moves so that ones moving towards a surface are more likely to be accepted."""

    def __init__(self, weight):
        self.weight = weight
        super(TargetEdge, self).__init__(weight=weight)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )


class AvoidEdge(Bias):
    """Bias next-moves so that ones away from a surface are more likely to be accepted."""

    def __init__(self, weight):
        self.weight = weight
        super(AvoidEdge, self).__init__(weight=weight)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )


class TargetPath(Bias):
    """Bias next-moves so that ones following a pre-defined path are more likely to be accepted."""

    def __init__(self, target_path, weight):
        self.target_path = target_path
        self.weight = weight
        super(TargetPath, self).__init__(weight=weight)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )
