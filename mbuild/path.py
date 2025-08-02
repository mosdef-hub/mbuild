"""Molecular paths and templates"""

import math
from abc import abstractmethod

import freud
import numpy as np
from scipy.interpolate import interp1d

from mbuild import Compound
from mbuild.utils.geometry import bounding_box

"""
A path is basically a bond graph with coordinates/positions
assigned to the nodes. This is kind of what Compound is already.

The interesting and challenging part is building up/creating the path.
This follows an algorithm to generate next coordinates.
Any random path generation algorithm will include
a rejection/acception step. We basically end up with
Monte Carlo. Some path algorithms won't be random (lamellae)

Is Path essentially going to be a simple Monte carlo-ish
class that others can inherit from then implement their own approach?

Classes that inherit from path will have their own
verions of next_coordinate(), check_path(), etc..
We can define abstract methods for these in Path.
We can put universally useful methods in Path as well (e.g., n_list(), to_compound(), etc..).

Some paths (lamellar, cyclic, spiral, etc..) would kind of just do
everything in generate() without having to use
next_coordinate() or check_path(). These would still need to be
defined, but just left empty and/or always return True in the case of check_path.
Maybe that means these kinds of "paths" need a different data structure?

Do we have RandomPath and DeterministicPath?

RandomPath ideas:
- Random walk (tons of possibilities here)
- Branching
- Multiple self-avoiding random walks
- ??

DeterministicPath ideas:
- Lamellar layers
- Cyclic polymers
- Helix
- Spiral
- Knots

Some combination of these?
- Lamellar + random walk to generate semi-crystalline like structures?
- Make a new path by adding together multiple paths
- Some kind of data structure/functionality for new_path = Path(start_from_path=other_path)

Other TODOs:
    Make coordinates a property with a setter? Keep as a plain attribute?
"""

try:
    from numba import njit

    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False


class Path:
    def __init__(self, N=None, coordinates=None, bond_graph=None):
        self.bond_graph = bond_graph
        # Only N is defined, make empty coordinates array with size N
        # Use case: Random walks
        if N is not None and coordinates is None:
            self.N = N
            self.coordinates = np.zeros((N, 3))
        # Only coordinates is defined, set N from length
        # Use case: class method from_coordinates
        elif coordinates is not None and N is None:
            self.N = len(coordinates)
            self.coordinates = coordinates
        # Neither is defined, use list for coordinates
        # Use case: Lamellar - Won't know N initially
        elif N is None and coordinates is None:
            self.N = N
            self.coordinates = []
        else:
            raise ValueError("Specify either one of N and coordinates, or neither")
        self.generate()
        if self.N is None:
            self.N = len(self.coordinates)

    @classmethod
    def from_coordinates(cls, coordinates, bond_graph=None):
        return cls(coordinates=coordinates, bond_graph=bond_graph, N=None)

    def _extend_coordinates(self, N):
        zeros = np.zeros((N, 3))
        new_array = np.concatenate(self.coordinates, zeros)
        self.coordinates = new_array
        self.N += N

    @abstractmethod
    def generate(self):
        """Abstract class for running a Path generation algorithm

        This method should:
        -----------------
            - Set initial conditions
            - Implement Path generation steps by calling _next_coordinate() and _check_path()
            - Handle cases of next coordiante acceptance
            - Handle cases of next coordinate rejection
        """
        pass

    @abstractmethod
    def _next_coordinate(self):
        """Algorithm to generate the next coordinate in the path"""
        pass

    @abstractmethod
    def _check_path(self):
        """Algorithm to accept/reject trial move of the current path"""
        pass

    def neighbor_list(self, r_max, coordinates=None, box=None):
        """Use freud to create a neighbor list of a set of coordinates."""
        if coordinates is None:
            coordinates = self.coordinates
        if box is None:
            box = bounding_box(coordinates)
        freud_box = freud.box.Box(Lx=box[0], Ly=box[1], Lz=box[2])
        aq = freud.locality.AABBQuery(freud_box, coordinates)
        aq_query = aq.query(
            query_points=coordinates,
            query_args=dict(r_min=0.0, r_max=r_max, exclude_ii=True),
        )
        nlist = aq_query.toNeighborList()
        return nlist

    def to_compound(self, bead_name="_A", bead_mass=1):
        """Visualize a path as an mBuild Compound"""
        compound = Compound()
        for xyz in self.coordinates:
            compound.add(Compound(name=bead_name, mass=bead_mass, pos=xyz))
        if self.bond_graph:
            compound.set_bond_graph(self.bond_graph)
        return compound

    def apply_mapping(self):
        # TODO: Finish, add logic to align orientation with path site pos and bond graph
        """Mapping other compounds onto a Path's coordinates

        mapping = {"A": "c1ccccc1C=C", "B": "C=CC=C"}
        """
        pass

    def _path_history(self):
        """Maybe this is a method that can be used optionally.
        We could add a save_history parameter to __init__.
        Depending on the approach, saving histories might add additional
        computation time and resources.
        Might be useful for more complicated random walks/branching algorithms
        """
        pass


class HardSphereRandomWalk(Path):
    def __init__(
        self,
        N,
        bond_length,
        radius,
        min_angle,
        max_angle,
        max_attempts,
        seed,
        trial_batch_size=5,
        start_from_path=None,
        start_from_path_index=None,
        tolerance=1e-5,
        use_numba=False,
        bond_graph=None,
    ):
        """Generates coordinates from a self avoiding random walk using
        fixed bond lengths, hard spheres, and minimum and maximum angles
        formed by 3 consecutive points.

        Possible angle values are sampled uniformly between min_angle
        and max_angle between the new site and the two previous sites.

        Parameters:
        -----------
        bond_length : float, required
            Fixed bond length between 2 coordinates.
        radius : float, required
            Radius of sites used in checking for overlaps.
        min_angle : float, required
            Minimum angle used when randomly selecting angle
            for the next step.
        max_angle : float, required
            Maximum angle used when randomly selecting angle
            for the next step.
        seed : int, default = 42
            Random seed
        tolerance : float, default = 1e-4
            Tolerance used for rounding.
        bond_graph : networkx.graph.Graph; optional
            Sets the bonding of sites along the path.
        """
        self.bond_length = bond_length
        self.radius = radius
        self.min_angle = min_angle
        self.max_angle = max_angle
        self.seed = seed
        self.tolerance = tolerance
        self.trial_batch_size = int(trial_batch_size)
        self.max_attempts = int(max_attempts)
        self.attempts = 0
        self.start_from_path_index = start_from_path_index
        self.start_from_path = start_from_path
        if start_from_path and start_from_path_index is not None:
            # TODO: Do we need np.copy here?
            coordinates = np.concatenate(
                (
                    np.copy(start_from_path.coordinates).astype(np.float32),
                    np.zeros((N, 3), dtype=np.float32),
                ),
                axis=0,
            )
            self.count = len(start_from_path.coordinates) - 1
            N = None
        else:
            coordinates = np.zeros((N, 3), dtype=np.float32)
            self.count = 0
            self.start_index = 0
        # Get methods to use for random walk
        if not use_numba:
            self.next_coordinate = _random_coordinate_numpy
            self.check_path = _check_path_numpy
        elif use_numba and _NUMBA_AVAILABLE:
            self.next_coordinate = _random_coordinate_numba
            self.check_path = _check_path_numba
        elif use_numba and not _NUMBA_AVAILABLE:
            raise RuntimeError(
                "numba was not found. Set `use_numba` to `False` or add numba to your environment."
            )

        super(HardSphereRandomWalk, self).__init__(
            coordinates=coordinates, N=None, bond_graph=bond_graph
        )

    def generate(self):
        np.random.seed(self.seed)
        if not self.start_from_path:
            # With fixed bond lengths, first move is always accepted
            phi = np.random.uniform(0, 2 * np.pi)
            theta = np.random.uniform(0, np.pi)
            next_pos = np.array(
                [
                    self.bond_length * np.sin(theta) * np.cos(phi),
                    self.bond_length * np.sin(theta) * np.sin(phi),
                    self.bond_length * np.cos(theta),
                ]
            )
            self.coordinates[1] = next_pos
            self.count += 1  # We already have 1 accepted move
        else:
            # Start a while loop here
            started_next_path = False
            while not started_next_path:
                new_xyzs = self.next_coordinate(
                    pos1=self.start_from_path.coordinates[self.start_from_path_index],
                    pos2=self.start_from_path.coordinates[
                        self.start_from_path_index - 1
                    ],
                    bond_length=self.bond_length,
                    min_angle=self.min_angle,
                    max_angle=self.max_angle,
                    batch_size=self.trial_batch_size,
                )
                new_xyz_found = False
                for xyz in new_xyzs:
                    if self.check_path(
                        existing_points=self.coordinates[: self.count + 1],
                        new_point=xyz,
                        radius=self.radius,
                        tolerance=self.tolerance,
                    ):
                        self.coordinates[self.count + 1] = xyz
                        self.count += 1
                        self.attempts += 1
                        new_xyz_found = True
                        started_next_path = True
                        break
                if not new_xyz_found:
                    self.attempts += 1

                if self.attempts == self.max_attempts and self.count < self.N:
                    raise RuntimeError(
                        "The maximum number attempts allowed have passed, and only ",
                        f"{self.count} sucsessful attempts were completed.",
                        "Try changing the parameters and running again.",
                    )

        while self.count < self.N - 1:
            new_xyzs = self.next_coordinate(
                pos1=self.coordinates[self.count],
                pos2=self.coordinates[self.count - 1],
                bond_length=self.bond_length,
                min_angle=self.min_angle,
                max_angle=self.max_angle,
                batch_size=self.trial_batch_size,
            )
            for xyz in new_xyzs:
                if self.check_path(
                    existing_points=self.coordinates[: self.count + 1],
                    new_point=xyz,
                    radius=self.radius,
                    tolerance=self.tolerance,
                ):
                    self.coordinates[self.count + 1] = xyz
                    self.count += 1
                    break
            self.attempts += 1

            if self.attempts == self.max_attempts and self.count < self.N:
                raise RuntimeError(
                    "The maximum number attempts allowed have passed, and only ",
                    f"{self.count} sucsessful attempts were completed.",
                    "Try changing the parameters and running again.",
                )


class Lamellar(Path):
    def __init__(
        self,
        num_layers,
        layer_separation,
        layer_length,
        bond_length,
        num_stacks=1,
        stack_separation=None,
        bond_graph=None,
    ):
        self.num_layers = num_layers
        self.layer_separation = layer_separation
        self.layer_length = layer_length
        self.bond_length = bond_length
        self.num_stacks = num_stacks
        self.stack_separation = stack_separation
        super(Lamellar, self).__init__(N=None, bond_graph=bond_graph)

    def generate(self):
        layer_spacing = np.arange(0, self.layer_length, self.bond_length)
        # Info needed for generating coords of the curves between layers
        r = self.layer_separation / 2
        arc_length = r * np.pi
        arc_num_points = math.floor(arc_length / self.bond_length)
        arc_angle = np.pi / (arc_num_points + 1)  # incremental angle
        arc_angles = np.linspace(arc_angle, np.pi, arc_num_points, endpoint=False)
        #        stack_coordinates = []
        for i in range(self.num_layers):
            if i % 2 == 0:  # Even layer; build from left to right
                layer = [
                    np.array([self.layer_separation * i, y, 0]) for y in layer_spacing
                ]
                # Mid-point between this and next layer; use to get curve coords.
                origin = layer[-1] + np.array([r, 0, 0])
                arc = [
                    origin + np.array([-np.cos(theta), np.sin(theta), 0]) * r
                    for theta in arc_angles
                ]
            else:  # Odd layer; build from right to left
                layer = [
                    np.array([self.layer_separation * i, y, 0])
                    for y in layer_spacing[::-1]
                ]
                # Mid-point between this and next layer; use to get curve coords.
                origin = layer[-1] + np.array([r, 0, 0])
                arc = [
                    origin + np.array([-np.cos(theta), -np.sin(theta), 0]) * r
                    for theta in arc_angles
                ]
            if i != self.num_layers - 1:
                self.coordinates.extend(layer + arc)
            else:
                self.coordinates.extend(layer)
        if self.num_stacks > 1:
            first_stack_coordinates = np.copy(np.array(self.coordinates))
            # Now find info for curves between stacked layers
            r = self.stack_separation / 2
            arc_length = r * np.pi
            arc_num_points = math.floor(arc_length / self.bond_length)
            arc_angle = np.pi / (arc_num_points + 1)  # incremental angle
            arc_angles = np.linspace(arc_angle, np.pi, arc_num_points, endpoint=False)
            if self.num_layers % 2 == 0:
                odd_stack_mult = -1
            else:
                odd_stack_mult = 1
            for i in range(1, self.num_stacks):
                if i % 2 != 0:  # Odd stack
                    this_stack = np.copy(first_stack_coordinates[::-1]) + np.array(
                        [0, 0, self.stack_separation * i]
                    )
                    origin = self.coordinates[-1] + np.array([0, 0, r])
                    arc = [
                        origin
                        + np.array([0, odd_stack_mult * np.sin(theta), np.cos(theta)])
                        * r
                        for theta in arc_angles
                    ]
                    self.coordinates.extend(arc[::-1])
                    self.coordinates.extend(list(this_stack))
                elif i % 2 == 0:  # Even stack
                    this_stack = np.copy(first_stack_coordinates) + np.array(
                        [0, 0, self.stack_separation * i]
                    )
                    origin = self.coordinates[-1] + np.array([0, 0, r])
                    arc = [
                        origin + np.array([0, -np.sin(theta), np.cos(theta)]) * r
                        for theta in arc_angles
                    ]
                    self.coordinates.extend(arc[::-1])
                    self.coordinates.extend(list(this_stack))

    def _next_coordinate(self):
        pass

    def _check_path(self):
        pass


class StraightLine(Path):
    def __init__(self, spacing, N, direction=(1, 0, 0), bond_graph=None):
        self.spacing = spacing
        self.direction = np.asarray(direction)
        super(StraightLine, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        self.coordinates = np.array(
            [np.zeros(3) + i * self.spacing * self.direction for i in range(self.N)]
        )

    def _next_coordinate(self):
        pass

    def _check_path(self):
        pass


class Cyclic(Path):
    def __init__(self, spacing=None, N=None, radius=None, bond_graph=None):
        self.spacing = spacing
        self.radius = radius
        n_params = sum(1 for i in (spacing, N, radius) if i is not None)
        if n_params != 2:
            raise ValueError("You must specify only 2 of spacing, N and radius.")
        super(Cyclic, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        if self.spacing and self.N:
            self.radius = (self.N * self.spacing) / (2 * np.pi)
        elif self.radius and self.spacing:
            self.N = int((2 * np.pi * self.radius) / self.spacing)
        else:
            self.spacing = (2 * np.pi) / self.N

        angles = np.arange(0, 2 * np.pi, (2 * np.pi) / self.N)
        self.coordinates = np.array(
            [(np.cos(a) * self.radius, np.sin(a) * self.radius, 0) for a in angles]
        )

    def _next_coordinate(self):
        pass

    def _check_path(self):
        pass


class Knot(Path):
    def __init__(self, spacing, N, m, bond_graph=None):
        self.spacing = spacing
        self.m = m
        super(Knot, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        t_dense = np.linspace(0, 2 * np.pi, 5000)
        # Base (unscaled) curve
        if self.m == 3:  # Trefoil knot (3_1)
            R, r = 1.0, 0.3
            x = (R + r * np.cos(3 * t_dense)) * np.cos(2 * t_dense)
            y = (R + r * np.cos(3 * t_dense)) * np.sin(2 * t_dense)
            z = r * np.sin(3 * t_dense)
        elif self.m == 4:  # Figure-eight knot (4_1)
            x = (2 + np.cos(2 * t_dense)) * np.cos(3 * t_dense)
            y = (2 + np.cos(2 * t_dense)) * np.sin(3 * t_dense)
            z = np.sin(4 * t_dense)
        elif self.m == 5:  # Cinquefoil knot (5_1), a (5,2) torus knot
            R, r = 1.0, 0.3
            x = (R + r * np.cos(5 * t_dense)) * np.cos(2 * t_dense)
            y = (R + r * np.cos(5 * t_dense)) * np.sin(2 * t_dense)
            z = r * np.sin(5 * t_dense)
        else:
            raise ValueError("Only m=3, m=4 and m=5 are supported.")
        # Compute arc length of a base curve
        coords_dense = np.stack((x, y, z), axis=1)
        deltas = np.diff(coords_dense, axis=0)
        dists = np.linalg.norm(deltas, axis=1)
        arc_lengths = np.concatenate([[0], np.cumsum(dists)])
        base_length = arc_lengths[-1]
        L_target = (self.N - 1) * self.spacing
        # Scale to match target contour length
        scale = L_target / base_length
        coords_dense *= scale
        arc_lengths *= scale
        # Resample uniformly along arc length based on target separation and N sites
        desired_arcs = np.linspace(0, L_target, self.N, endpoint=False)
        x_interp = interp1d(arc_lengths, coords_dense[:, 0])(desired_arcs)
        y_interp = interp1d(arc_lengths, coords_dense[:, 1])(desired_arcs)
        z_interp = interp1d(arc_lengths, coords_dense[:, 2])(desired_arcs)
        self.coordinates = np.stack((x_interp, y_interp, z_interp), axis=1)


class Helix(Path):
    def __init__(
        self, N, radius, rise, twist, right_handed=True, bottom_up=True, bond_graph=None
    ):
        """
        Generate helical path.

        Parameters:
        -----------
        radius : float, required
            Radius of the helix (nm)
        rise : float, required
            Rise per site on path (nm)
        twist : float, required
            Twist per site in path (degrees)
        right_handed : bool, default True
            Set the handedness of the helical twist
            Set to false for a left handed twist
        bottom_up : bool, default True
            If True, the twist is in the positive Z direction
            If False, the twist is in the negative Z direction

        Notes:
        ------
        To create a double helix pair (e.g., DNA) create two paths
        with opposite values for right_handed and bottom_up
        """
        self.radius = radius
        self.rise = rise
        self.twist = twist
        self.right_handed = right_handed
        self.bottom_up = bottom_up
        super(Helix, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        indices = reversed(range(self.N)) if not self.bottom_up else range(self.N)
        for i in indices:
            angle = np.deg2rad(i * self.twist)
            if not self.right_handed:
                angle *= -1
            x = self.radius * np.cos(angle)
            y = self.radius * np.sin(angle)
            z = i * self.rise if self.bottom_up else -i * self.rise
            self.coordinates[i] = (x, y, z)


class Spiral2D(Path):
    def __init__(self, N, a, b, spacing, bond_graph=None):
        """
        Generate a 2D spiral path in the XY plane.

        Parameters
        ----------
        N: int, required
            Number of sites in the path
        a : float, required
            The initial radius (nm)
        b : float, required
            Determines how fast radius grows per angle increment (r = a + bθ)
        spacing : float, required
            Distance between adjacent sites (nm)
        """
        self.a = a
        self.b = b
        self.spacing = spacing
        super().__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        theta = 0.0
        for i in range(self.N):
            r = self.a + self.b * theta
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            z = 0.0
            self.coordinates[i] = (x, y, z)
            # Estimate next angle increment based on arc length
            ds_dtheta = np.sqrt((r) ** 2 + self.b**2)
            dtheta = self.spacing / ds_dtheta
            theta += dtheta


class ZigZag(Path):
    def __init__(
        self,
        N,
        spacing=1.0,
        angle_deg=120.0,
        sites_per_segment=1,
        plane="xy",
        bond_graph=None,
    ):
        self.spacing = spacing
        self.angle_deg = angle_deg
        self.sites_per_segment = sites_per_segment
        self.plane = plane
        if N % sites_per_segment != 0:
            raise ValueError("N must be evenly divisible by sites_per_segment")
        super(ZigZag, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        angle_rad = np.deg2rad(self.angle_deg)
        direction = np.array([1.0, 0.0])
        position = np.zeros(2)
        coords = []
        step_count = 0
        segment_count = 0

        for i in range(self.N):
            coords.append(position.copy())
            position += self.spacing * direction
            step_count += 1

            # Rotate
            if step_count == self.sites_per_segment:
                step_count = 0
                sign = -1 if segment_count % 2 == 0 else 1
                rot_matrix = np.array(
                    [
                        [np.cos(sign * angle_rad), -np.sin(sign * angle_rad)],
                        [np.sin(sign * angle_rad), np.cos(sign * angle_rad)],
                    ]
                )
                direction = rot_matrix @ direction
                segment_count += 1

        # Map into 3D space based on chosen plane
        for i, (x2d, y2d) in enumerate(coords):
            if self.plane == "xy":
                self.coordinates[i] = (x2d, y2d, 0)
            elif self.plane == "xz":
                self.coordinates[i] = (x2d, 0, y2d)
            elif self.plane == "yz":
                self.coordinates[i] = (0, x2d, y2d)


# Internal helper/utility methods below:


def _random_coordinate_numpy(pos1, pos2, bond_length, min_angle, max_angle, batch_size):
    # Vector formed by previous 2 coordinates
    v1 = pos2 - pos1
    v1_norm = v1 / np.linalg.norm(v1)
    # Generate batch of random angles
    thetas = np.random.uniform(min_angle, max_angle, size=batch_size)
    # Batch of random vectors and center around origin (0,0,0)
    r = np.random.rand(batch_size, 3) - 0.5
    dot_products = np.dot(r, v1_norm)
    r_perp = r - dot_products[:, np.newaxis] * v1_norm
    norms = np.linalg.norm(r_perp, axis=1)
    for norm in norms:
        if norm < 1e-6:
            norm = 1.0
    r_perp_norm = r_perp / norms[:, np.newaxis]
    v2s = (
        np.cos(thetas)[:, np.newaxis] * v1_norm
        + np.sin(thetas)[:, np.newaxis] * r_perp_norm
    )
    next_positions = pos1 + v2s * bond_length
    return next_positions


def _check_path_numpy(existing_points, new_point, radius, tolerance):
    """Check new trial point against previous ones only."""
    sq_dists = np.sum((existing_points - new_point) ** 2, axis=1)
    min_sq_dist = (radius - tolerance) ** 2
    return not np.any(sq_dists < min_sq_dist)


@njit(fastmath=True)
def norm(vec):
    s = 0.0
    for i in range(vec.shape[0]):
        s += vec[i] * vec[i]
    return np.sqrt(s)


@njit(fastmath=True)
def _random_coordinate_numba(
    pos1,
    pos2,
    bond_length,
    min_angle,
    max_angle,
    batch_size,
):
    v1 = pos2 - pos1
    v1_norm = v1 / norm(v1)

    thetas = np.empty(batch_size, dtype=np.float32)
    for i in range(batch_size):
        thetas[i] = min_angle + (max_angle - min_angle) * np.random.random()

    r = np.empty((batch_size, 3), dtype=np.float32)
    for i in range(batch_size):
        for j in range(3):
            r[i, j] = np.random.random() - 0.5

    dot_products = np.empty(batch_size, dtype=np.float32)
    for i in range(batch_size):
        dot = 0.0
        for j in range(3):
            dot += r[i, j] * v1_norm[j]
        dot_products[i] = dot

    r_perp = np.empty((batch_size, 3), dtype=np.float32)
    for i in range(batch_size):
        for j in range(3):
            r_perp[i, j] = r[i, j] - dot_products[i] * v1_norm[j]

    norms = np.empty(batch_size, dtype=np.float32)
    for i in range(batch_size):
        norms[i] = norm(r_perp[i])

    for i in range(batch_size):
        if norms[i] < 1e-6:
            norms[i] = 1.0

    r_perp_norm = np.empty((batch_size, 3), dtype=np.float32)
    for i in range(batch_size):
        for j in range(3):
            r_perp_norm[i, j] = r_perp[i, j] / norms[i]

    v2s = np.empty((batch_size, 3), dtype=np.float32)
    for i in range(batch_size):
        cos_theta = np.cos(thetas[i])
        sin_theta = np.sin(thetas[i])
        for j in range(3):
            v2s[i, j] = cos_theta * v1_norm[j] + sin_theta * r_perp_norm[i, j]

    next_positions = np.empty((batch_size, 3), dtype=np.float32)
    for i in range(batch_size):
        for j in range(3):
            next_positions[i, j] = pos1[j] + v2s[i, j] * bond_length
    return next_positions


@njit
def _check_path_numba(existing_points, new_point, radius, tolerance):
    min_sq_dist = (radius - tolerance) ** 2
    for i in range(existing_points.shape[0]):
        dist_sq = 0.0
        for j in range(existing_points.shape[1]):
            diff = existing_points[i, j] - new_point[j]
            dist_sq += diff * diff
        if dist_sq < min_sq_dist:
            return False
    return True
