"""Classes to generate intra-molecular paths and configurations."""

import math
from abc import abstractmethod
from copy import deepcopy

import networkx as nx
import numpy as np
from scipy.interpolate import interp1d

from mbuild import Compound
from mbuild.path.path_utils import check_path, random_coordinate


class Path:
    """Creates a path from a given set of coordinates and a bond graph.
    This class is designed to be use in mbuild.polymer.Polymer.build_from_path().
    This also serves as the base class from which other paths inherit from.

    Parameters
    ----------
    N : int, optional
        The number of sites belonging to the path.
    coordinates : array-like, optional
        Creates a path from a pre-defined set of coordinates
    bond_graph : networkx.Graph, optional
       The graph defining the edges between coordinates
    """

    def __init__(self, N=None, coordinates=None, bond_graph=None, bead_name="_A"):
        self.bond_graph = bond_graph
        self.bead_name = bead_name
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
        # Use case: Lamellar - Don't know N initially
        elif N is None and coordinates is None:
            self.coordinates = []
        else:
            raise ValueError("Specify either one of N and coordinates, or neither.")
        self.__post_init__()

    def __post_init__(self):
        """Needed for CodeQL in order to call abstract method inside __init__()"""
        self.generate()
        try:
            self.N
        except:
            self.N = len(self.coordinates)
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

    def add_edge(self, u, v):
        """Add an edge to the Path's bond graph.
        This sets attributes for the edge which include: bond direction, bond length and bond type
        u -> v corresponds to previous site and current site where the bond direction is calculated as v - u.
        """
        bond_vec = self.coordinates[v] - self.coordinates[u]
        bond_length = np.linalg.norm(bond_vec)
        bond_vec /= bond_length
        # Get node names from previous step to current step
        u_name = self.bond_graph.nodes[u]["name"]
        v_name = self.bond_graph.nodes[v]["name"]
        self.bond_graph.add_edge(
            u_of_edge=u,
            v_of_edge=v,
            direction=bond_vec.tolist(),
            length=float(bond_length),
            bond_type=(u_name, v_name),
        )

    def get_bonded_sites(self):
        """Get all bonded pairs and their bond-vector orientations."""
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )

    def get_coordinates(self):
        if isinstance(self.coordinates, list):
            return np.array(self.coordinates)
        else:
            return self.coordinates

    @abstractmethod
    def generate(self):
        """Abstract class for running a Path generation algorithm.
        Sub-classes that inherit from Path should implement their
        own method under genrate.
        """
        pass

    def to_compound(self):
        """Convert a path and its bond graph to an mBuild Compound."""
        compound = Compound()
        for node_id, attrs in self.bond_graph.nodes(data=True):
            compound.add(Compound(name=attrs["name"], pos=attrs["xyz"]))
        compound.set_bond_graph(self.bond_graph)
        return compound

    def apply_mapping(self):
        # TODO: Finish, add logic to align orientation with path site pos and bond graph
        """Mapping other compounds onto a Path's coordinates

        mapping = {"A": "c1ccccc1C=C", "B": "C=CC=C"}
        mapping = {"A": mb.Compound, "B": mb.Compound}

        for bond in bond graph edges:
            site u: add compound
            rotate site u head-tail vec to align with bond direction
            set orientation and separation for site u head port
            rotate site v tail-head vec to align with bond direction
            set orientation and separation for site v tail port

        """
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )

    def _path_history(self):
        """Maybe this is a method that can be used optionally.
        We could add a save_history parameter to __init__.
        Depending on the approach, saving histories might add additional
        computation time and resources.
        Might be useful for more complicated random walks/branching algorithms
        """
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )


class HardSphereRandomWalk(Path):
    def __init__(
        self,
        N,
        bond_length,
        radius,
        min_angle,
        max_angle,
        bead_name="_A",
        volume_constraint=None,
        start_from_path=None,
        start_from_path_index=None,
        attach_paths=False,
        initial_point=None,
        include_compound=None,
        max_attempts=1e5,
        seed=42,
        trial_batch_size=20,
        tolerance=1e-5,
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
            Minimum angle (radians) used when randomly selecting angle
            for the next step.
        max_angle : float, required
            Maximum angle (radians) used when randomly selecting angle
            for the next step.
        bead_name : str, default = "_A"
            The name assigned to each site in this random walk.
        volume_constraint : mbuild.utils.volumes.Constraint, optional
            Used to reject moves which are outside of the volume constraint
        start_from_path : mbuild.path.Path, optional
            An instance of a previous Path to start the random walk from.
            This path's sites are used in checking for overlapping sites.
        start_from_path_index : int, optional
            An index of `start_from_path` used as the starting point for this random walk.
        attach_paths : bool, default = False
            If True, adds an edge between the starting point of the last path and the first
            point of this path.
        initial_point : array-like, optional
            Used as the coordinate for the first site in this random walk path.
        include_compound : mb.Compound, optional
            An mBuild compound that is considered in this random walk.
            This random walk with reject new points that overlap with the Compound's coordinates.
        seed : int, default = 42
            Random seed
        trial_batch_size : int, default = 5
            The number of trial moves to attempt in parallel for each step.
            Using larger values can improve success rates for more dense
            random walks.
        max_attempts : int, default = 1e5
            The maximum number of trial moves to attempt before quiting.
            for the random walk.
        tolerance : float, default = 1e-4
            Tolerance used for rounding and checking for overlaps.

        Notes
        -----
        Each next move can be attempted in batches, set by the ``trial_batch_size``
        parameter. The batch size moves do not count towards the maximum allowed
        attempts. For example, 1 random walk step with a trail batch size of 20 counts at
        only one attempted move. Larger values of ``trial_batch_size`` may help
        highly constrained walks finish, but may hurt performance.

        You can start a random walk from a previously created path with the
        ``start_from_path`` and ``start_from_path_index`` parameters. For example,
        create a ``Lamellar`` path and run a random walk from its last index to begin
        generating a semi-crystalline like structure. Or, string together multiple
        HardSphereRandomWalk paths where the final result of each is passed into the
        next random walk.
        """
        if initial_point is not None:
            self.initial_point = np.asarray(initial_point)
        else:
            self.initial_point = None
        self.include_compound = include_compound
        self.bond_length = bond_length
        self.radius = radius
        self.min_angle = min_angle
        self.max_angle = max_angle
        self.seed = seed
        self.bead_name = bead_name
        self.volume_constraint = volume_constraint
        self.tolerance = tolerance
        self.trial_batch_size = int(trial_batch_size)
        self.max_attempts = int(max_attempts)
        self.attempts = 0
        self.start_from_path_index = start_from_path_index
        self.start_from_path = start_from_path
        self.attach_paths = attach_paths
        self._particle_pairs = {}

        # This random walk is including a previous path
        if start_from_path:
            coordinates = np.concatenate(
                (
                    start_from_path.get_coordinates().astype(np.float32),
                    np.zeros((N, 3), dtype=np.float32),
                ),
                axis=0,
            )
            self.count = len(start_from_path.coordinates)
            N = None
            bond_graph = deepcopy(start_from_path.bond_graph)
            if start_from_path_index is not None and start_from_path_index < 0:
                self.start_from_path_index = self.count + start_from_path_index
        else:  # Not starting from another path
            bond_graph = nx.Graph()
            coordinates = np.zeros((N, 3), dtype=np.float32)
            N = None
            self.count = 0
            self.start_index = 0
        # Need this for error message about reaching max tries
        self._init_count = self.count
        # Select methods to use for random walk
        # Hard-coded for now, possible to make other RW methods and pass them in
        self.next_step = random_coordinate
        self.check_path = check_path
        # Create RNG state.
        self.rng = np.random.default_rng(seed)
        super(HardSphereRandomWalk, self).__init__(
            coordinates=coordinates, N=N, bond_graph=bond_graph
        )

    def generate(self):
        initial_xyz = self._initial_points()
        # Set the first coordinate
        self.coordinates[self.count] = initial_xyz
        self.bond_graph.add_node(
            self.count,
            name=self.bead_name,
            xyz=self.coordinates[self.count],
        )
        if not self.start_from_path and not self.volume_constraint:
            # If no volume constraint, then first move is always accepted
            phi = self.rng.uniform(0, 2 * np.pi)
            theta = self.rng.uniform(0, np.pi)
            next_pos = np.array(
                [
                    self.bond_length * np.sin(theta) * np.cos(phi),
                    self.bond_length * np.sin(theta) * np.sin(phi),
                    self.bond_length * np.cos(theta),
                ]
            )
            self.coordinates[1] = self.coordinates[0] + next_pos
            self.count += 1
            self.bond_graph.add_node(
                self.count,
                name=self.bead_name,
                xyz=self.coordinates[self.count],
            )
            self.add_edge(u=self.count - 1, v=self.count)
        # Not starting from another path, but have a volume constraint
        # Possible for second point to be out-of-bounds
        elif not self.start_from_path and self.volume_constraint:
            self.coordinates[0] = initial_xyz
            next_point_found = False
            while not next_point_found:
                phi = self.rng.uniform(0, 2 * np.pi)
                theta = self.rng.uniform(0, np.pi)
                xyz = np.array(
                    [
                        self.bond_length * np.sin(theta) * np.cos(phi),
                        self.bond_length * np.sin(theta) * np.sin(phi),
                        self.bond_length * np.cos(theta),
                    ]
                )
                is_inside_mask = self.volume_constraint.is_inside(
                    points=np.array([xyz]), buffer=self.radius
                )
                if np.all(is_inside_mask):
                    self.coordinates[1] = self.coordinates[0] + xyz
                    self.count += 1
                    self.bond_graph.add_node(
                        self.count,
                        name=self.bead_name,
                        xyz=self.coordinates[self.count],
                    )
                    self.add_edge(u=self.count - 1, v=self.count)
                    self.attempts += 1
                    next_point_found = True
                # 2nd point failed, continue while loop
                self.attempts += 1

                if self.attempts == self.max_attempts and self.count < self.N:
                    raise RuntimeError(
                        "The maximum number attempts allowed have passed, and only ",
                        f"{self.count - self._init_count} sucsessful attempts were completed.",
                        "Try changing the parameters or seed and running again.",
                    )

        # Starting random walk from a previous set of coordinates (another path)
        # This point was accepted in self._initial_point with these conditions
        # If attach_paths, then add edge between first node of this path and node of last path
        else:
            self.bond_graph.add_node(self.count, name=self.bead_name, xyz=initial_xyz)
            if self.attach_paths:
                self.add_edge(u=self.start_from_path_index, v=self.count)
            self.attempts += 1

        # Initial conditions set (points 1 and 2), now start RW with min/max angles
        while self.count < self.N - 1:
            # Choosing angles and vectors is the only random part
            # Get a batch of these once, pass them into numba functions
            batch_angles, batch_vectors = self._generate_random_trials()
            new_xyzs = self.next_step(
                pos1=self.coordinates[self.count],
                pos2=self.coordinates[self.count - 1],
                bond_length=self.bond_length,
                thetas=batch_angles,
                r_vectors=batch_vectors,
                batch_size=self.trial_batch_size,
            )
            if self.volume_constraint:
                is_inside_mask = self.volume_constraint.is_inside(
                    points=new_xyzs, buffer=self.radius
                )
                new_xyzs = new_xyzs[is_inside_mask]
            for xyz in new_xyzs:
                if self.include_compound:
                    # Include compound particle coordinates in check for overlaps
                    existing_points = np.concatenate(
                        (self.coordinates[: self.count + 1], self.include_compound.xyz)
                    )
                else:
                    existing_points = self.coordinates[: self.count + 1]
                # Now check for overlaps for each trial point
                # Stop after the first success
                if self.check_path(
                    existing_points=existing_points,
                    new_point=xyz,
                    radius=self.radius,
                    tolerance=self.tolerance,
                ):
                    self.coordinates[self.count + 1] = xyz
                    self.count += 1
                    self.bond_graph.add_node(self.count, name=self.bead_name, xyz=xyz)
                    self.add_edge(u=self.count - 1, v=self.count)
                    break
            self.attempts += 1

            if self.attempts == self.max_attempts and self.count < self.N:
                raise RuntimeError(
                    "The maximum number attempts allowed have passed, and only ",
                    f"{self.count - self._init_count} sucsessful attempts were completed.",
                    "Try changing the parameters or seed and running again.",
                )

    def _generate_random_trials(self):
        """Generate a batch of random angles and vectors using the RNG state."""
        # Batch of angles used to determine next positions
        thetas = self.rng.uniform(
            self.min_angle, self.max_angle, size=self.trial_batch_size
        ).astype(np.float32)
        # Batch of random vectors and center around origin (0,0,0)
        r = self.rng.uniform(-0.5, 0.5, size=(self.trial_batch_size, 3)).astype(
            np.float32
        )
        return thetas, r

    def _initial_points(self):
        """Choosing first and second points depends on the parameters passed in."""
        # Use manually specified initial point
        if self.initial_point is not None:
            return self.initial_point

        # Random initial point, no volume constraint: Bounds set by radius and N steps
        elif not any(
            [self.volume_constraint, self.initial_point, self.start_from_path]
        ):
            max_dist = (self.N * self.radius) - self.radius
            xyz = self.rng.uniform(low=-max_dist, high=max_dist, size=3)
            return xyz

        # Random point inside volume constraint, not starting from another path
        elif self.volume_constraint and not self.start_from_path_index:
            xyz = self.rng.uniform(
                low=self.volume_constraint.mins + self.radius,
                high=self.volume_constraint.maxs - self.radius,
                size=3,
            )
            return xyz

        # Starting from another path, run Monte Carlo
        # Accepted next move is first point of this random walk
        elif self.start_from_path and self.start_from_path_index is not None:
            # TODO: handle start_from_path index of negative values
            # Set to the corresponding actual index value of the last path
            if self.start_from_path_index == 0:
                pos2_coord = 1  # Use the second (1) point of the last path for angles
            else:  # use the site previous to start_from_path_index for angles
                pos2_coord = self.start_from_path_index - 1
            started_next_path = False
            while not started_next_path:
                batch_angles, batch_vectors = self._generate_random_trials()
                new_xyzs = self.next_step(
                    pos1=self.start_from_path.get_coordinates()[
                        self.start_from_path_index
                    ],
                    pos2=self.start_from_path.get_coordinates()[pos2_coord],
                    bond_length=self.bond_length,
                    thetas=batch_angles,
                    r_vectors=batch_vectors,
                    batch_size=self.trial_batch_size,
                )
                if self.volume_constraint:
                    is_inside_mask = self.volume_constraint.is_inside(
                        points=new_xyzs, buffer=self.radius
                    )
                    new_xyzs = new_xyzs[is_inside_mask]
                for xyz in new_xyzs:
                    if self.include_compound:
                        existing_points = np.concatenate(
                            (
                                self.coordinates[: self.count + 1],
                                self.include_compound.xyz,
                            )
                        )
                    else:
                        existing_points = self.coordinates[: self.count + 1]
                    if self.check_path(
                        existing_points=existing_points,
                        new_point=xyz,
                        radius=self.radius,
                        tolerance=self.tolerance,
                    ):
                        return xyz
                self.attempts += 1

                if self.attempts == self.max_attempts and self.count < self.N:
                    raise RuntimeError(
                        "The maximum number attempts allowed have passed, and only ",
                        f"{self.count - self._init_count} sucsessful attempts were completed.",
                        "Try changing the parameters or seed and running again.",
                    )


class Lamellar(Path):
    """Generate a 2-D or 3-D lamellar-like path.

    Parameters
    ----------
    spacing : float (nm), required
        The distance between two adjacent sites in the path.
    num_layers : int, required
        The number of times the lamellar path curves around
        creating another layer.
    layer_separation : float (nm), required
        The distance between any two layers.
    num_stacks : int, required
        The number of times to repeat each layer in the Z direction.
        Setting this to 1 creates a single, 2D lamellar-like path.
    stack_separation : float (nm), required
        The distance between two stacked layers.
    """

    def __init__(
        self,
        num_layers,
        layer_separation,
        layer_length,
        bond_length,
        num_stacks=1,
        bead_name="_A",
        stack_separation=None,
    ):
        self.num_layers = num_layers
        self.layer_separation = layer_separation
        self.layer_length = layer_length
        self.bond_length = bond_length
        self.num_stacks = num_stacks
        self.stack_separation = stack_separation
        bond_graph = nx.Graph()
        super(Lamellar, self).__init__(
            N=None, bond_graph=bond_graph, bead_name=bead_name
        )

    def generate(self):
        layer_spacing = np.arange(
            0, self.layer_length, self.bond_length, dtype=np.float64
        )
        # Info needed for generating coords of the arc curves between layers
        r = self.layer_separation / 2
        arc_length = r * np.pi
        arc_num_points = math.floor(arc_length / self.bond_length)
        arc_angle = np.pi / (arc_num_points + 1)  # incremental angle
        arc_angles = np.linspace(arc_angle, np.pi, arc_num_points, endpoint=False)
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
            else:  # Last layer, don't include another arc set of coordinates
                self.coordinates.extend(layer)
        if self.num_stacks > 1:
            first_stack_coordinates = np.copy(np.array(self.coordinates))
            # Get info for curves between stacked layers
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
        # Create linear (path) bond graph
        for i, xyz in enumerate(self.coordinates):
            self.bond_graph.add_node(
                i,
                name=self.bead_name,
                xyz=xyz,
            )
            if i != 0:
                self.add_edge(u=i - 1, v=i)


class StraightLine(Path):
    """Generates a set of coordinates in a straight line along a given axis.

    Parameters
    ----------
    spacing : float, required
        The distance between sites along the path.
    N : int, required
        The number of sites in the path.
    direction : array-like (1,3), default = (1,0,0)
        The direction to align the straight path along.
    """

    def __init__(self, spacing, N, direction=(1, 0, 0), bond_graph=None):
        self.spacing = spacing
        self.direction = np.asarray(direction)
        super(StraightLine, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        self.coordinates = np.array(
            [np.zeros(3) + i * self.spacing * self.direction for i in range(self.N)]
        )


class Cyclic(Path):
    """Generates a set of coordinates evenly spaced along a circle.

    Parameters
    ----------
    spacing : float, optional
        Distance between sites along the path.
    N : int, optional
        Number of sites in the cyclic path.
    radius : float, optional
        The radius (nm) of the cyclic path.

    Notes
    -----
    Only two of spacing, N and radius can be defined, as the third
    is determined by the other two.

    If using this Path to build a cyclic polymer, be sure to
    set ``bond_head_tail = True`` in ``mbuild.polymer.Polymer.build_from_path``
    """

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


class Knot(Path):
    """Generate a knot path.

    Parameters
    ----------
    spacing : float (nm)
        The spacing between sites along the path.
    N : int
        The number of total sites in the path.
    m : int in [3, 4, 5]
        The number of crossings in the knot.
        3 gives the trefoil knot, 4 gives the figure 8 knot and 5 gives the cinquefoil knot.
        Only values of 3, 4 and 5 are currently supported.
    bond_graph : networkx.graph
        Sets the bond graph between sites.
    """

    def __init__(self, spacing, N, m, bond_graph=None):
        self.spacing = spacing
        self.m = m
        super(Knot, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        # Generate dense sites first, sample actual ones later from spacing
        # Prevents spacing between sites changing with curvature
        t_dense = np.linspace(0, 2 * np.pi, 5000)
        # Base (unscaled) curve
        if self.m == 3:  # Trefoil knot  https://en.wikipedia.org/wiki/Trefoil_knot
            R, r = 1.0, 0.3
            x = (R + r * np.cos(3 * t_dense)) * np.cos(2 * t_dense)
            y = (R + r * np.cos(3 * t_dense)) * np.sin(2 * t_dense)
            z = r * np.sin(3 * t_dense)
        elif (
            self.m == 4
        ):  # Figure-eight https://en.wikipedia.org/wiki/Figure-eight_knot_(mathematics)
            x = (2 + np.cos(2 * t_dense)) * np.cos(3 * t_dense)
            y = (2 + np.cos(2 * t_dense)) * np.sin(3 * t_dense)
            z = np.sin(4 * t_dense)
        elif (
            self.m == 5
        ):  # Cinquefoil knot https://en.wikipedia.org/wiki/Cinquefoil_knot
            R, r = 1.0, 0.3
            x = (R + r * np.cos(5 * t_dense)) * np.cos(2 * t_dense)
            y = (R + r * np.cos(5 * t_dense)) * np.sin(2 * t_dense)
            z = r * np.sin(5 * t_dense)
        else:
            raise ValueError("Only m=3, m=4 and m=5 are currently supported.")
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
        """Generate helical path.

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
            Set to False for a left handed twist
        bottom_up : bool, default True
            If True, the twist is in the positive Z direction
            If False, the twist is in the negative Z direction

        Notes:
        ------
        To create a double helix pair (e.g., DNA) create two paths
        with opposite values for right_handed and bottom_up.
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
        """Generate a 2D spiral path in the XY plane.

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
    """Generates a path following a zig-zag pattern in a given plane.

    Parameters
    ----------
    spacing : float, default = 1.0
        The distance between consecutive sites along the path.
    angle_deg : float, default = 120.
        The rotation applied between segments
    sites_per_segment : int, default = 4
        The number of sites before rotating and beginning next segment.
    plane : str, default = "xy"
        The plane that the sites in the path occupy
    bond_graph : networkx.Graph
        Defines connectivity between sites

    """

    def __init__(
        self,
        N,
        spacing=1.0,
        angle_deg=120.0,
        sites_per_segment=4,
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
