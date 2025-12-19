import time

import numpy as np


class Termination:
    def __init__(self, criteria):
        self.criteria = list(criteria)
        self.required = [c for c in self.criteria if c.required_to_end]
        self.not_required = [c for c in self.criteria if not c.required_to_end]
        # TODO, keep a list of triggered critera, add to logging when walk ends
        self.triggered = []

    def attach_path(self, path):
        for c in self.criteria:
            c._attach_path(path)

    def is_met(self):
        # Check required criteria first
        if all([c.is_met() for c in self.required]):
            return True
        if any([c.is_met() for c in self.not_required]):
            return True
        return False


class Criterion:
    """"""

    def __init__(self, required_to_end):
        self.required_to_end = required_to_end

    def _attach_path(self, path):
        self.path = path
        self.rng = getattr(path, "rng", None)

    def is_met(self):
        """Implemented in sub classes of Criterion."""
        raise NotImplementedError


class NumSites(Criterion):
    def __init__(self, num_sites, required_to_end=True):
        self.num_sites = num_sites
        super().__init__(required_to_end)

    def is_met(self):
        return self.path.count - self.path.init_count >= self.num_sites


class NumAttempts(Criterion):
    def __init__(self, max_attempts, required_to_end=False):
        self.max_attempts = int(max_attempts)
        super().__init__(required_to_end)

    def is_met(self):
        return self.path.attempts >= self.max_attempts


class WallTime(Criterion):
    def __init__(self, max_time, required_to_end=False):
        self.max_time = max_time
        super().__init__(required_to_end)

    def is_met(self):
        current_time = time.time()
        total_time = current_time - self.path.start_time
        return total_time >= self.max_time


class WithinCoordinate(Criterion):
    def __init__(
        self, final_coordinate, distance, tolerance=1e-3, required_to_end=True
    ):
        self.distance = float(distance)
        self.final_coordinate = np.asarray(final_coordinate)
        self.tolerance = tolerance
        super().__init__(required_to_end)

    def is_met(self):
        last_site = self.path.coordinates[self.path.count]
        current_distance = np.linalg.norm(self.final_coordinate - last_site)
        return current_distance <= self.distance + self.tolerance


class FinalCoordinate(Criterion):
    def __init__(self, final_coordinate, tolerance=1e-3, required_to_end=True):
        self.final_coordinate = np.asarray(final_coordinate)
        self.tolerance = tolerance
        super().__init__(required_to_end)

    def is_met(self):
        last_site = self.path.coordinates[self.path.count]
        return np.linalg.norm(self.final_coordinate - last_site) <= self.tolerance


class PairDistance(Criterion):
    def __init__(self, distance, pair_type=None, required_to_end=True):
        self.distance = float(distance)
        self.pair_type = pair_type
        super().__init__(required_to_end)

    def is_met(self):
        raise NotImplementedError


class RadiusGyration(Criterion):
    def __init__(self, radius_of_gyration, tolerance=1e-2, required_to_end=True):
        self.radius_of_gyration = float(radius_of_gyration)
        self.tolerance = tolerance
        self.max_rg2 = (radius_of_gyration + tolerance) ** 2
        super().__init__(required_to_end)

    def is_met(self):
        # Enforce at least 3 sites
        n = self.path.count - self.path.init_count + 1
        if n < 2:
            return False

        coords = self.path.coordinates[self.path.init_count : self.path.count + 1]
        center = coords.mean(axis=0)
        diffs = coords - center
        rg2 = np.mean(np.sum(diffs * diffs, axis=1))
        return rg2 >= self.max_rg2


class EndtoEndDistance(Criterion):
    def __init__(self, distance, tolerance=1e-2, required_to_end=True):
        self.distance = float(distance)
        self.tolerance = tolerance
        super().__init__(required_to_end)

    def is_met(self):
        # Enforce at least 3 sites
        n = self.path.count - self.path.init_count + 1
        if n < 2:
            return False
        first_site = self.path.coordinates[self.path.init_count]
        last_site = self.path.coordinates[self.path.count]
        return np.linalg.norm(last_site - first_site) >= self.distance - self.tolerance
