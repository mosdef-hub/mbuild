import time

import numpy as np


class Termination:
    def __init__(self, terminators):
        if isinstance(terminators, Terminator):
            self.terminators = [terminators]
        else:
            self.terminators = list(terminators)
        # These must all be True to trigger termination
        self.required_to_end = [i for i in self.terminators if i.required_to_end]
        # Don't need to be True, but are used as safe-guards (WallTime, NumAttempts)
        self.not_required_to_end = [
            i for i in self.terminators if not i.required_to_end
        ]
        # TODO, keep a list of triggered critera, add to logging when walk ends
        self.triggered = []
        self.success = False

    def _attach_path(self, path):
        """This is automatically called within HardSphereRandomWalk."""
        for i in self.terminators:
            i._attach_path(path)

    def is_met(self):
        # Check required criteria first
        if all([i.is_met() for i in self.required_to_end]):
            self.success = True
            return True
        if any([i.is_met() for i in self.not_required_to_end]):
            return True
        return False

    def summarize(self):
        """Print a quick summary of termination status."""
        lines = []
        status = "SUCCESS" if self.success else "ABORTED"
        lines.append(f"Termination status: {status}")
        lines.append("")

        for term in self.terminators:
            name = term.__class__.__name__
            met = term.is_met()
            role = "required" if term.required_to_end else "safeguard"
            flag = "✓" if met else "✗"
            lines.append(f"[{flag}] {name} ({role})")
        return "\n".join(lines)


class Terminator:
    """"""

    def __init__(self, required_to_end):
        self.required_to_end = required_to_end

    def _attach_path(self, path):
        self.path = path
        self.rng = getattr(path, "rng", None)

    def is_met(self):
        """Implemented in sub classes of Terminator."""
        raise NotImplementedError


class NumSites(Terminator):
    def __init__(self, num_sites, required_to_end=True):
        self.num_sites = num_sites
        super().__init__(required_to_end)

    def is_met(self):
        return self.path.count - self.path._init_count >= self.num_sites - 1


class NumAttempts(Terminator):
    def __init__(self, max_attempts, required_to_end=False):
        self.max_attempts = int(max_attempts)
        super().__init__(required_to_end)

    def is_met(self):
        return self.path.attempts >= self.max_attempts


class WallTime(Terminator):
    def __init__(self, max_time, required_to_end=False):
        self.max_time = max_time
        super().__init__(required_to_end)

    def is_met(self):
        current_time = time.time()
        total_time = current_time - self.path.start_time
        return total_time >= self.max_time


class WithinCoordinate(Terminator):
    def __init__(
        self, target_coordinate, distance, tolerance=1e-3, required_to_end=True
    ):
        self.distance = float(distance)
        self.target_coordinate = np.asarray(target_coordinate)
        self.tolerance = tolerance
        super().__init__(required_to_end)

    def is_met(self):
        last_site = self.path.coordinates[self.path.count]
        current_distance = np.linalg.norm(self.target_coordinate - last_site)
        return current_distance <= self.distance + self.tolerance


class PairDistance(Terminator):
    def __init__(self, distance, pair_type=None, required_to_end=True):
        self.distance = float(distance)
        self.pair_type = pair_type
        super().__init__(required_to_end)

    def is_met(self):
        raise NotImplementedError


class RadiusOfGyration(Terminator):
    def __init__(self, radius_of_gyration, tolerance=0.01, required_to_end=True):
        self.radius_of_gyration = float(radius_of_gyration)
        self.tolerance = tolerance
        self.rg2 = (radius_of_gyration + tolerance) ** 2
        super().__init__(required_to_end)

    def is_met(self):
        # Enforce at least 3 sites
        n = self.path.count - self.path._init_count + 1
        if n < 2:
            return False
        coords = self.path.coordinates[self.path._init_count : self.path.count + 1]
        center = coords.mean(axis=0)
        diffs = coords - center
        rg2 = np.mean(np.sum(diffs * diffs, axis=1))
        return self.rg2 - self.tolerance <= rg2 <= self.rg2 + self.tolerance


class EndToEndDistance(Terminator):
    def __init__(self, distance, tolerance=0.01, required_to_end=True):
        self.distance = float(distance)
        self.tolerance = tolerance
        super().__init__(required_to_end)

    def is_met(self):
        # Enforce at least 3 sites
        n = self.path.count - self.path._init_count + 1
        if n < 2:
            return False
        first_site = self.path.coordinates[self.path._init_count]
        last_site = self.path.coordinates[self.path.count]
        return (
            self.distance - self.tolerance
            <= np.linalg.norm(last_site - first_site)
            <= self.distance + self.tolerance
        )
