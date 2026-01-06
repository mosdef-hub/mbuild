"""Modular and composable rules for termination an mbuild.path.HardSphereRandomWalk."""

import time

import numpy as np


class Termination:
    """A modular and composable container for individual Terminator instances.
    This class should be passed to mbuild.path.HardSphereRandomWalk to provide
    instructions on when to end a random walk.

    Parameters
    ----------
    terminators : list, required
        A list-like object of mbuild.path.Terminator subclasses.
        Each are checked after each step of a random walk.

    Methods
    -------
    is_met(): Called internally in HardSphereRandomWalk.
        Returns `True` in two cases:

        1) Every `Terminator` instance with `is_target=True` is met.
        2) Any `Terminator` instance with `is_target=False` is met.

    summarize(): Prints a summary of all `Terminator` instances and thier status

    """

    def __init__(self, terminators):
        if isinstance(terminators, Terminator):
            self.terminators = [terminators]
        else:
            self.terminators = list(terminators)
        # These must all be True to trigger termination
        self.is_target = [i for i in self.terminators if i.is_target]
        # Don't need to be True, but are used as safe-guards (WallTime, NumAttempts)
        self.not_is_target = [i for i in self.terminators if not i.is_target]
        # TODO, keep a list of triggered critera, add to logging when walk ends
        self.triggered = []
        self.success = False

    def _attach_path(self, path):
        """This is automatically called within HardSphereRandomWalk."""
        for i in self.terminators:
            i._attach_path(path)

    def is_met(self):
        # Check required criteria first
        if all([i.is_met() for i in self.is_target]):
            self.success = True
            return True
        if any([i.is_met() for i in self.not_is_target]):
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
            role = "required" if term.is_target else "safeguard"
            flag = "✓" if met else "✗"
            lines.append(f"[{flag}] {name} ({role})")
        return "\n".join(lines)


class Terminator:
    """Defines a single condition that can trigger a termination of a HardSphereRandomWalk.

    Parameters
    ----------
    is_target : bool, required
        If `True`, this terminator represents a target condition that the walk
        is attempting to reach. Triggering a target condition indicates
        successful completion of the walk. If `False`, the terminator represents
        a safeguard that limits execution (e.g., maximum attempts or wall time);
        triggering a safeguard causes the walk to stop without being considered
        successful.

    Notes
    -----
    Multiple instances of `Terminator` can be passed to `Termination`. Target
    terminators define completion goals (e.g., placing a fixed number of sites),
    while non-target terminators act as safety limits when convergence is
    unlikely (e.g., maximum attempts or total wall time).

    """

    def __init__(self, is_target):
        self.is_target = is_target

    def _attach_path(self, path):
        self.path = path
        self.rng = getattr(path, "rng", None)

    def is_met(self):
        """Implemented in sub classes of Terminator."""
        raise NotImplementedError


class NumSites(Terminator):
    """A terminator that triggers after a certain number of successful steps.

    Parameters
    ----------
    num_sites : int, required
        The number of successful sites before triggering.
    is_target : bool, default True
        If `True`, this terminator represents a target condition that the walk
        is attempting to reach. Triggering a target condition indicates
        successful completion of the walk. If `False`, the terminator represents
        a safeguard that limits execution (e.g., maximum attempts or wall time);
        triggering a safeguard causes the walk to stop without being considered
        successful.

    """

    def __init__(self, num_sites, is_target=True):
        self.num_sites = num_sites
        super().__init__(is_target)

    def is_met(self):
        return self.path.count - self.path._init_count >= self.num_sites - 1


class NumAttempts(Terminator):
    """A terminator that triggers after a certain number of total attempts.

    Parameters
    ----------
    max_attempts : int, required
        The number of attempts (successful or not) before triggering.
    is_target : bool, default False
        If `True`, this terminator represents a target condition that the walk
        is attempting to reach. Triggering a target condition indicates
        successful completion of the walk. If `False`, the terminator represents
        a safeguard that limits execution (e.g., maximum attempts or wall time);
        triggering a safeguard causes the walk to stop without being considered
        successful.

    """

    def __init__(self, max_attempts, is_target=False):
        self.max_attempts = int(max_attempts)
        super().__init__(is_target)

    def is_met(self):
        return self.path.attempts >= self.max_attempts


class WallTime(Terminator):
    """A terminator that triggers after a certain total time (seconds) as ellapsed.

    Parameters
    ----------
    max_time : int, required
        The total amount of time (seconds) allowed to ellapse before triggering.
    is_target : bool, default False
        If `True`, this terminator represents a target condition that the walk
        is attempting to reach. Triggering a target condition indicates
        successful completion of the walk. If `False`, the terminator represents
        a safeguard that limits execution (e.g., maximum attempts or wall time);
        triggering a safeguard causes the walk to stop without being considered
        successful.

    """

    def __init__(self, max_time, is_target=False):
        self.max_time = max_time
        super().__init__(is_target)

    def is_met(self):
        current_time = time.time()
        total_time = current_time - self.path.start_time
        return total_time >= self.max_time


class WithinCoordinate(Terminator):
    """A terminator that triggers after a successful step is within a cutoff distance to a target coordinate.

    Parameters
    ----------
    target_coordinate : array-like (n,3), required
        The target coordinate in units of nm.
    distance : float, required
        A cutoff distance (nm) from the target coordinate.
    is_target : bool, default True
        If `True`, this terminator represents a target condition that the walk
        is attempting to reach. Triggering a target condition indicates
        successful completion of the walk. If `False`, the terminator represents
        a safeguard that limits execution (e.g., maximum attempts or wall time);
        triggering a safeguard causes the walk to stop without being considered
        successful.

    """

    def __init__(self, target_coordinate, distance, tolerance=1e-3, is_target=True):
        self.distance = float(distance)
        self.target_coordinate = np.asarray(target_coordinate)
        self.tolerance = tolerance
        super().__init__(is_target)

    def is_met(self):
        last_site = self.path.coordinates[self.path.count]
        current_distance = np.linalg.norm(self.target_coordinate - last_site)
        return current_distance <= self.distance + self.tolerance


class PairDistance(Terminator):
    """A terminator that triggers after a successful step is within a cutoff distance to a certain site type.

    Parameters
    ----------
    pair_type : str or list of str, required
        The name of the site type(s) to consider in distance calculations
    distance : float, required
        A cutoff distance (nm) from the target sites.
    is_target : bool, default True
        If `True`, this terminator represents a target condition that the walk
        is attempting to reach. Triggering a target condition indicates
        successful completion of the walk. If `False`, the terminator represents
        a safeguard that limits execution (e.g., maximum attempts or wall time);
        triggering a safeguard causes the walk to stop without being considered
        successful.

    """

    def __init__(self, distance, pair_type=None, is_target=True):
        self.distance = float(distance)
        self.pair_type = pair_type
        super().__init__(is_target)

    def is_met(self):
        raise NotImplementedError


class RadiusOfGyration(Terminator):
    """A terminator that triggers after reaching a target square radius of gyration.

    Parameters
    ----------
    radius_of_gyration : float, required
        The target square radius of gyration (Rg^2) in units of nm^2.
    tolerance : float, default 0.01
        The allowable tolerance relative to the target. In units of nm^2.
    is_target : bool, default True
        If `True`, this terminator represents a target condition that the walk
        is attempting to reach. Triggering a target condition indicates
        successful completion of the walk. If `False`, the terminator represents
        a safeguard that limits execution (e.g., maximum attempts or wall time);
        triggering a safeguard causes the walk to stop without being considered
        successful.

    """

    def __init__(self, radius_of_gyration, tolerance=0.01, is_target=True):
        self.rg = float(radius_of_gyration)
        self.tolerance = tolerance
        super().__init__(is_target)

    def is_met(self):
        # Enforce at least 3 sites
        n = self.path.count - self.path._init_count + 1
        if n < 2:
            return False
        coords = self.path.coordinates[self.path._init_count : self.path.count + 1]
        center = coords.mean(axis=0)
        diffs = coords - center
        rg2 = np.mean(np.sum(diffs * diffs, axis=1))
        return self.rg - self.tolerance <= rg2 <= self.rg + self.tolerance


class EndToEndDistance(Terminator):
    """A terminator that triggers after reaching a target end-to-end distance.

    Parameters
    ----------
    distance : float, required
        The target end-to-end distance in units of nm.
    tolerance : float, default 0.01
        The allowable tolerance relative to the target. In units of nm.
    is_target : bool, default True
        If `True`, this terminator represents a target condition that the walk
        is attempting to reach. Triggering a target condition indicates
        successful completion of the walk. If `False`, the terminator represents
        a safeguard that limits execution (e.g., maximum attempts or wall time);
        triggering a safeguard causes the walk to stop without being considered
        successful.

    """

    def __init__(self, distance, tolerance=0.01, is_target=True):
        self.distance = float(distance)
        self.tolerance = tolerance
        super().__init__(is_target)

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
