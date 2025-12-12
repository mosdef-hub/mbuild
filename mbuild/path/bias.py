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

    def _attach_path(self, path):
        self.path = path
        # Inherit rng from the path for use in Bias classes
        self.rng = self.path.rng

    def __call__(self, candidates):
        raise NotImplementedError


class TargetCoordinate(Bias):
    """Bias next-moves so that ones moving closer to a final coordinate are more likely to be accepted."""

    def __init__(self, target_coordinate, weight):
        self.target_coordinate = np.asarray(target_coordinate)
        super(TargetCoordinate, self).__init__(weight=weight)

    def __call__(self, candidates):
        sq_distances = target_sq_distances(self.target_coordinate, candidates)
        # Large beta diminishes the effect of noise
        beta = self.weight / max(1e-6, (1.0 - self.weight))
        noise_scale = 1 - self.weight
        noise = self.rng.normal(0, noise_scale, size=sq_distances.shape)
        scores = beta * sq_distances + noise
        # Target coordinate should favor short distances, sort in ascending order (np default)
        sort_idx = np.argsort(scores)
        return candidates[sort_idx]


class AvoidCoordinate(Bias):
    """Bias next-moves so that ones moving further from a specific coordinate are more likely to be accepted."""

    def __init__(self, avoid_coordinate, weight):
        self.avoid_coordinate = np.asarray(avoid_coordinate)
        super(AvoidCoordinate, self).__init__(weight=weight)

    def __call__(self, candidates):
        sq_distances = target_sq_distances(
            self.avoid_coordinate.astype(np.float32), candidates.astype(np.float32)
        )
        # Large beta diminishes the effect of noise
        beta = self.weight / max(1e-6, (1.0 - self.weight))
        noise_scale = 1 - self.weight
        noise = self.rng.normal(0, noise_scale, size=sq_distances.shape)
        # Use posiive beta term to favor larger distances
        scores = beta * sq_distances + noise
        # Avoid cooardinate should favor larger distances, sort in descending order
        sort_idx = np.argsort(scores)[::-1]
        return candidates[sort_idx]


class TargetType(Bias):
    """Bias next-moves so that ones moving towards a specific site type are more likely to be accepted."""

    def __init__(self, target_type, weight, r_cut):
        self.target_type = target_type
        self.r_cut = r_cut
        super(TargetType, self).__init__(weight=weight)

    def __call__(self, candidates):
        types = np.array(
            [node[1]["name"] for node in self.path.bond_graph.nodes(data=True)]
        )
        # path.coordinates is an array with place holder values for future sites
        # Only check as far as types are defined (i.e., once a node is added)
        target_coords = self.path.coordinates[: len(types)][
            types == self.target_type
        ].astype(np.float32)
        densities = target_density(
            candidates=candidates, target_coords=target_coords, r_cut=self.r_cut
        )
        # Apply weights
        beta = self.weight / max(1e-6, (1.0 - self.weight))
        noise_scale = 1 - self.weight
        noise = self.rng.normal(0, noise_scale, size=densities.shape)
        scores = beta * densities + noise
        # Target type should favor larger densities, sort in descending order
        sort_idx = np.argsort(scores)[::-1]
        return candidates[sort_idx]


class AvoidType(Bias):
    """Bias next-moves so that ones moving away from a specific site type are more likely to be accepted."""

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
        # Apply weights
        beta = self.weight / max(1e-6, (1.0 - self.weight))
        noise_scale = 1 - self.weight
        noise = self.rng.normal(0, noise_scale, size=densities.shape)
        scores = beta * densities + noise
        # Avoid type should facor smaller densities, sort in ascending order (np default)
        sort_idx = np.argsort(scores)
        return candidates[sort_idx]


class TargetDirection(Bias):
    """Bias next-moves so that ones moving along a certain direction are more likely to be accepted."""

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
        beta = self.weight / max(1e-6, (1.0 - self.weight))
        noise_scale = 1.0 - self.weight
        noise = self.rng.normal(0.0, noise_scale, size=alignment.shape)
        scores = beta * alignment + noise
        # Larger dot product = better alignment with target, sort descending
        sort_idx = np.argsort(scores)[::-1]
        return candidates[sort_idx]


class AvoidDirection(Bias):
    """Bias next-moves so that ones not moving along a certain direction are more likely to be accepted."""

    def __init__(self, direction, weight):
        self.direction = direction
        super(AvoidDirection, self).__init__()

    def __call__(self, candidates):
        last_step_pos = self.path.coordinates[self.path.count]
        next_step_vectors = candidates - last_step_pos
        norms = np.linalg.norm(next_step_vectors, axis=1, keepdims=True)
        # Avoid divisions by zero
        norms = np.clip(norms, 1e-12, None)
        next_step_unit_vectors = next_step_vectors / norms
        # Alignment score: dot product with target direction
        alignment = np.dot(next_step_unit_vectors, self.direction)
        beta = self.weight / max(1e-6, (1.0 - self.weight))
        noise_scale = 1.0 - self.weight
        noise = self.rng.normal(0.0, noise_scale, size=alignment.shape)
        scores = beta * alignment + noise
        # Larger dot product = better alignment with target, sort ascending
        sort_idx = np.argsort(scores)
        return candidates[sort_idx]


class TargetEdge(Bias):
    """Bias next-moves so that ones moving towards a surface are more likely to be accepted."""

    def __init__(self, weight, system_coordinates, volume_constraint, new_coordinates):
        self.volume_constraint = volume_constraint
        super(TargetEdge, self).__init__(system_coordinates, new_coordinates)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )


class AvoidEdge(Bias):
    """Bias next-moves so that ones away from a surface are more likely to be accepted."""

    def __init__(self, weight, system_coordinates, volume_constraint, new_coordinates):
        self.volume_constraint = volume_constraint
        super(AvoidEdge, self).__init__(system_coordinates, new_coordinates)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )


class TargetPath(Bias):
    """Bias next-moves so that ones following a pre-defined path are more likely to be accepted."""

    def __init__(self, target_path, weight, system_coordinates, new_coordinates):
        self.target_path = target_path
        super(TargetPath, self).__init__(system_coordinates, new_coordinates)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )
