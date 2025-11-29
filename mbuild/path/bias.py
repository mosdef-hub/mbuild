"""Contains biases that can be included in mbuild.path.HardSphereRandomWalk."""

import numpy as np

from mbuild.path.path_utils import target_sq_distances


class Bias:
    def __init__(self, new_coordinates, path_coordinates=None, site_types=None):
        # Batch of candidate coordinates generated in Path
        self.new_coordinates = new_coordinates
        # Complete system coords needed for density and direction biases
        self.path_coordinates = path_coordinates
        # List of site types corresponding to path_coordinates
        # Needed for TargetType and AvoidType
        self.site_types = site_types
        self.N = len(self.system_coordinates)

    def __call__(self):
        raise NotImplementedError


class TargetCoordinate(Bias):
    """Bias next-moves so that ones moving closer to a final coordinate are more likely to be accepted."""

    def __init__(self, target_coordinate, weight, new_coordinates):
        self.target_coordinate = target_coordinate
        super(TargetCoordinate, self).__init__(new_coordinates=new_coordinates)

    def __call__(self):
        sq_distances = target_sq_distances(self.target_coordinate, self.new_coordinates)
        sort_idx = np.argsort(sq_distances)
        return self.new_points[sort_idx]


class AvoidCoordinate(Bias):
    """Bias next-moves so that ones moving further from a specific coordinate are more likely to be accepted."""

    def __init__(self, avoid_coordinate, weight, new_coordinates):
        self.avoid_coordinate = avoid_coordinate
        super(AvoidCoordinate, self).__init__(new_coordinates=new_coordinates)

    def __call__(self):
        sq_distances = target_sq_distances(self.avoid_coordinate, self.new_coordinates)
        # Sort in descending order, largest to smallest
        sort_idx = np.argsort(sq_distances)[::-1]
        return self.new_points[sort_idx]


class TargetType(Bias):
    """Bias next-moves so that ones moving towards a specific site type are more likely to be accepted."""

    def __init__(self, target_type, weight, path_coordinates, new_coordinates):
        self.target_type = target_type
        super(TargetType, self).__init__(new_coordinates=new_coordinates, path_coordinates=path_coordinates)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )


class AvoidType(Bias):
    """Bias next-moves so that ones moving away from a specific site type are more likely to be accepted."""

    def __init__(self, avoid_type, weight, path_coordinates, new_coordinates):
        self.avoid_type = avoid_type
        super(AvoidType, self).__init__(new_coordinates=new_coordinates, path_coordinates=path_coordinates)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )


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


class TargetDirection(Bias):
    """Bias next-moves so that ones moving along a certain direction are more likely to be accepted."""

    def __init__(self, direction, weight, system_coordinates, new_coordinates):
        self.direction = direction
        super(TargetDirection, self).__init__(system_coordinates, new_coordinates)

    def __call__(self):
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )


class AvoidDirection(Bias):
    """Bias next-moves so that ones not moving along a certain direction are more likely to be accepted."""

    def __init__(self, direction, weight, system_coordinates, new_coordinates):
        self.direction = direction
        super(AvoidDirection, self).__init__(system_coordinates, new_coordinates)

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
