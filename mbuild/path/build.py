"""Classes to generate intra-molecular paths and configurations."""

import logging
import math
import time

import networkx as nx
import numpy as np
from scipy.interpolate import interp1d

from mbuild import Compound
from mbuild.exceptions import PathConvergenceError
from mbuild.path.constraints import CuboidConstraint, CylinderConstraint
from mbuild.path.path_utils import (
    calculate_sq_distances,
    check_path,
    random_coordinate,
)
from mbuild.path.points import (
    AnglesSampler,
    generate_trials,
    get_initial_point,
    get_second_point,
)
from mbuild.path.termination import NumSites, Termination, Terminator

logger = logging.getLogger(__name__)


class Path:
    """Creates a path from a given set of coordinates and a bond graph.
    This class is designed to be use in mbuild.polymer.Polymer.build_from_path().

    Parameters
    ----------
    coordinates : array-like, optional
        Creates a path from a pre-defined set of coordinates
    bond_graph : networkx.Graph, optional
       The graph defining the edges between coordinates
    bead_name : str, default '_A'
        The name assigned to each site. This is helpful when using
        multiple `Path` instances to build heterogeneous systems.

    """

    # TODO, allow a sequence of bead_names
    def __init__(self, coordinates=None, bond_graph=None, bead_name="_A"):
        if (
            coordinates is not None
            and bond_graph is not None
            and isinstance(bead_name, np.ndarray)
        ):
            assert len(coordinates) == len(bond_graph), (
                len(coordinates),
                len(bond_graph),
            )
            assert len(coordinates) == len(bead_name), (
                len(coordinates),
                len(bead_name),
            )
            self.bond_graph = bond_graph
            self.coordinates = coordinates
            self.beads = bead_name
        elif coordinates is not None and bond_graph is not None:
            assert len(coordinates) == len(bond_graph)
            self.bond_graph = bond_graph
            self.coordinates = coordinates
            self.beads = np.array([bead_name for _ in range(len(coordinates))])
        elif coordinates is None and bond_graph is not None:
            self.bond_graph = bond_graph
            self.coordinates = np.array(  # might not have a bond_graph xyz
                [node.get("xyz") for node in bond_graph.nodes(data=True)]
            )
            self.beads = np.array(
                [node.get("name") for node in bond_graph.nodes(data=True)]
            )
        elif coordinates is not None and bond_graph is None:
            self.coordinates = np.asarray(coordinates)
            self.bond_graph = nx.Graph()
            self.beads = [bead_name] * (len(self.coordinates))
            for idx in range(len((self.coordinates))):
                self.bond_graph.add_node(idx)
        else:
            self.coordinates = np.array([], dtype=np.float32)
            self.bond_graph = nx.Graph()
            self.beads = np.array([], dtype="U10")

    def __eq__(self, other):
        return (
            np.all(self.coordinates == other.coordinates)
            and np.all(self.beads == other.beads)
            and nx.is_isomorphic(self.bond_graph, other.bond_graph)
        )

    def __add__(self, other):
        coordinates = np.concat((self.coordinates, other.coordinates))
        beads = np.concat((self.beads, other.beads))
        bond_graph = nx.compose(
            self.bond_graph, other.bond_graph
        )  # TODO: Don't overwrite nodes in bg
        return Path(coordinates, bond_graph, beads)

    @classmethod
    def from_compound(cls, compound):
        coordinates = compound.xyz

        # Create the path with coordinates and bond graph
        path = cls(coordinates=coordinates, bead_name=compound.name)
        path.bond_graph = nx.Graph()

        # Ensure all nodes have xyz and name attributes
        partDict = {}
        for idx, part in enumerate(compound.particles()):
            # Node doesn't exist in bond graph, add it
            path.bond_graph.add_node(idx, name=part.name, xyz=coordinates[idx])
            partDict[part] = idx

        # add edges
        for part1, part2 in compound.bonds():
            path.bond_graph.add_edge(partDict[part1], partDict[part2])

        return path

    def append_coordinates(self, points, bead_name="_A"):
        """Create new coordinates for appending values to."""
        if hasattr(points, "__len__") and len(points) == 3 and points.ndim == 1:
            points = np.asarray([points])
        if not isinstance(points, np.ndarray):
            try:
                points = np.ndarray(points)
            except TypeError:
                raise ValueError(f"{points=} must be an array of shape (N,3)")
        if points.ndim == 1:
            points = np.array([points])  # make a 2d array

        if self.coordinates.size == 0:
            self.coordinates = points
            self.beads = [bead_name] * len(points)
            self._extend_bond_graph(bead_name)
            return
        self.coordinates = np.concatenate((self.coordinates, points))
        self.beads = np.concatenate((self.beads, [bead_name] * len(points)))
        self._extend_bond_graph(bead_name)  # TODO: This won't need a bead name

    def _extend_coordinates(self, N):
        """Create new coordinates for setting values."""
        if self.coordinates.size == 0:
            self.coordinates = np.zeros((N, 3), dtype=np.float32)
            return
        zeros = np.zeros((N, 3), dtype=self.coordinates.dtype)
        new_array = np.concatenate([self.coordinates, zeros])
        self.coordinates = new_array

    def _extend_bond_graph(self, bead_name):
        """Make sure all coordinates are properly added to the bondgraph."""
        if len(self.bond_graph) == len(self.coordinates):
            return
        index = len(self.bond_graph)
        for coord in self.coordinates[len(self.bond_graph) :]:
            self.bond_graph.add_node(index, name=bead_name, xyz=coord)
            index += 1

    def _extend_beads(self, bead_name):
        """Update bead names to extended coordinates."""
        diff = len(self.coordinates) - len(self.beads)
        if not diff:
            return
        self.beads = np.concatenate((self.beads, [bead_name] * diff))

    def _connect_edges(self, connectivity, indices=None, attach_index=-1):
        """Adds edges to self.bond_graph matching a given style `connectivity`."""
        if indices is None:
            indices = np.arange(0, len(self.coordinates))

        if connectivity == "disconnected":
            return
        elif connectivity == "linear":
            for idx1, idx in zip(indices, indices[1:]):
                self.add_edge(idx1, idx)
        elif connectivity == "link-linear":
            for idx1, idx in zip(indices, indices[1:]):
                self.add_edge(idx1, idx)
            if attach_index in self.bond_graph.nodes:
                # link to previous node
                self.add_edge(attach_index, indices[0])
        elif connectivity == "cycle":
            for idx1, idx in zip(indices, indices[1:]):
                self.add_edge(idx1, idx)
            self.add_edge(indices[0], indices[-1])
        else:
            raise ValueError(
                f"Argument {connectivity=} is incorrect. Pass one of `linear`, `link-linear`, `cycle`"
            )

    def form_linear_bond_graph(self, indices=None):
        """Iterates through the Path's bond_graph by indices to build a linearly connected graph.

        Parameters
        ----------
        indices : list, default None
            iterate through the bond_graph by indices and form edges between nodes.


        Notes
        -----
        This assumes site i is always bonded to site i + 1 and i - 1 (i.e., linear graph).
        """
        n_coords = len(self.coordinates)
        if indices is None:  # assume all coords
            indices = list(range(n_coords))
        if indices[0] not in self.bond_graph:
            self.bond_graph.add_node(
                indices[0],
                name=self.beads[indices[0]],
                xyz=self.coordinates[indices[0]],
            )

        for idx1, idx in zip(indices, indices[1:]):
            if idx >= n_coords:
                raise ValueError(
                    f"Index {idx=} is out of bounds for Path with {n_coords}"
                )
            if idx not in self.bond_graph:
                self.bond_graph.add_node(
                    idx, name=self.beads[idx], xyz=self.coordinates[idx]
                )
                self.bond_graph.add_node(
                    idx, name=self.beads[idx], xyz=self.coordinates[idx]
                )
            self.add_edge(idx1, idx)

    def add_edge(self, u, v):
        """Add an edge to the Path's bond graph."""
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

    def find_neighbors(
        self, u, min_bond_length, max_bond_length, excluded_bond_depth=0
    ):
        """Don't add a new particle."""
        pass
        # return candidate_neighbors

    def get_bonded_sites(self):
        """Get all bonded pairs and their bond-vector orientations."""
        raise NotImplementedError(
            "This feature of mBuild 2.0 has not been implemented yet."
        )

    def to_compound(self):
        """Convert a path and its bond graph to an mBuild Compound."""
        compound = Compound()
        compounds = []
        for node_id, attrs in self.bond_graph.nodes(data=True):
            compounds.append(
                Compound(name=attrs["name"], pos=self.coordinates[node_id])
            )
        compound.add(compounds)
        for edge1, edge2 in self.bond_graph.edges():
            compound.bond_graph.add_edge(compounds[edge1], compounds[edge2])
        return compound

    def to_mol2(self):
        from mbuild.path.formats import to_mol2 as _to_mol2

        return _to_mol2(self)

    def to_mol(self):
        """
        Convert mBuild Path to SDF/MOL format

        Parameters:
        -----------
        path : mbuild.Path
        mol_name : str
            Name of the molecule
        atom_types : dict
            Mapping of atom indices to atom type symbols (e.g., {0: 'A', 5: 'B'})
            If None, all atoms will be 'A'
        """
        from mbuild.path.formats import to_mol as _to_mol

        return _to_mol(self)

    def to_mol3000(self, G=None):
        """
        Convert mBuild Path to SDF/MOL V3000 format

        Parameters:
        -----------
        G : nx.Graph, default None
            Bondgraph to use for visualization.
        """
        from mbuild.path.formats import to_mol3000 as _to_mol3000

        return _to_mol3000(self)

    def visualize(self, radius, hide_periodic_bonds=False):
        from mbuild.utils.visualize import visualize_path

        return visualize_path(self, radius, hide_periodic_bonds)

    def relax(self, bead_radius, bond_length=None, steps=1000, seed=1, nthreads=1):
        """Perform a dpd simulation to relax the current path."""
        # from mbuild.simulation import hoomd_cap_displacement, hoomd_fire, ForcesHandler, HoomdSimulation
        from mbuild.simulation import energy_minimize_path

        energy_minimize_path(self, bead_radius, bond_length, steps, seed, nthreads)
        return


def lamellar(
    path,
    num_layers,
    layer_separation,
    layer_length,
    bond_length,
    initial_point=(0, 0, 0),
    num_stacks=1,
    stack_separation=None,
    left_to_right=True,
    bead_name="_A",
):
    """Generate a 2-D or 3-D lamellar-like path.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The Path object to populate with coordinates
    num_layers : int, required
        The number of times the lamellar path curves around creating another layer.
    layer_separation : float (nm), required
        The distance between any two layers.
    layer_length : float (nm), required
        The distance of a lamellar layer before curving to the next.
    bond_length : float (nm), required
        The distance between two adjacent sites in the path.
    initial_point : nd.array (1,3), default (0,0,0)
        The coordinate of the first site of the lamellar path.
    num_stacks : int, default 1
        The number of times to repeat each layer in the Z direction.
    stack_separation : float (nm), optional
        The distance between two stacked layers. Required if `num_stacks` >= 2.
    left_to_right : boolean, default True
        If `True`, the first layer is built with increasing y-coordinates from the origin.
    """
    initial_point = np.asarray(initial_point)

    # Coordinates in the y-direction (layer-length) of the lamellar layer
    layer_spacing = np.arange(0, layer_length, bond_length)
    if not left_to_right:
        layer_spacing *= -1
    layer_spacing += initial_point[1]

    # Info needed for generating coords of the arc curves between layers
    r = layer_separation / 2
    arc_length = r * np.pi
    arc_num_points = math.floor(arc_length / bond_length)
    arc_angle = np.pi / (arc_num_points + 1)
    arc_angles = np.linspace(arc_angle, np.pi, arc_num_points, endpoint=False)

    coordinates = []

    # Iterate and build up layers
    for i in range(num_layers):
        x = initial_point[0] + (layer_separation * i)
        if i % 2 == 0:  # Even layer
            layer = [np.array([x, y, initial_point[2]]) for y in layer_spacing]
            origin = layer[-1] + np.array([r, 0, 0])
            if left_to_right:
                arc = [
                    origin + np.array([-np.cos(theta), np.sin(theta), 0]) * r
                    for theta in arc_angles
                ]
            else:
                arc = [
                    origin + np.array([-np.cos(theta), -np.sin(theta), 0]) * r
                    for theta in arc_angles
                ]
        else:  # Odd layer
            layer = [np.array([x, y, initial_point[2]]) for y in layer_spacing[::-1]]
            origin = layer[-1] + np.array([r, 0, 0])
            if left_to_right:
                arc = [
                    origin + np.array([-np.cos(theta), -np.sin(theta), 0]) * r
                    for theta in arc_angles
                ]
            else:
                arc = [
                    origin + np.array([-np.cos(theta), np.sin(theta), 0]) * r
                    for theta in arc_angles
                ]

        if i != num_layers - 1:
            coordinates.extend(layer + arc)
        else:
            coordinates.extend(layer)

    # Build up lamellar structure in 3rd dimension (Z) by stacking layers
    if num_stacks > 1:
        first_stack_coordinates = np.copy(np.array(coordinates))
        r = stack_separation / 2
        arc_length = r * np.pi
        arc_num_points = math.floor(arc_length / bond_length)
        arc_angle = np.pi / (arc_num_points + 1)
        arc_angles = np.linspace(arc_angle, np.pi, arc_num_points, endpoint=False)

        if num_layers % 2 == 0:
            if left_to_right:
                odd_stack_mult = -1
                even_stack_mult = -1
            else:
                odd_stack_mult = 1
                even_stack_mult = -1
        else:
            if left_to_right:
                odd_stack_mult = 1
                even_stack_mult = -1
            else:
                odd_stack_mult = -1
                even_stack_mult = 1

        for i in range(1, num_stacks):
            if i % 2 != 0:  # Odd stack
                this_stack = np.copy(first_stack_coordinates[::-1]) + np.array(
                    [0, 0, stack_separation * i]
                )
                origin = coordinates[-1] + np.array([0, 0, r])
                arc = [
                    origin
                    + np.array([0, odd_stack_mult * np.sin(theta), np.cos(theta)]) * r
                    for theta in arc_angles
                ]
                coordinates.extend(arc[::-1])
                coordinates.extend(list(this_stack))
            elif i % 2 == 0:  # Even stack
                this_stack = np.copy(first_stack_coordinates) + np.array(
                    [0, 0, stack_separation * i]
                )
                origin = coordinates[-1] + np.array([0, 0, r])
                arc = [
                    origin
                    + np.array([0, even_stack_mult * np.sin(theta), np.cos(theta)]) * r
                    for theta in arc_angles
                ]
                coordinates.extend(arc[::-1])
                coordinates.extend(list(this_stack))

    path.coordinates = np.array(coordinates)
    path.beads = [bead_name] * len(coordinates)
    path.form_linear_bond_graph()


def straight_line(path, spacing, N, direction=(1, 0, 0), bead_name="_A"):
    """Generates a set of coordinates in a straight line along a given axis.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The Path object to populate with coordinates
    spacing : float, required
        The distance between sites along the path.
    N : int, required
        The number of sites in the path.
    direction : array-like (1,3), default = (1,0,0)
        The direction to align the straight path along.
    """
    direction = np.asarray(direction)
    path.coordinates = np.array(
        [np.zeros(3) + i * spacing * direction for i in range(N)]
    )
    path.beads = [bead_name] * N
    path.form_linear_bond_graph()


def cyclic(path, spacing=None, N=None, radius=None, closed=True, bead_name="_A"):
    """Generates a set of coordinates evenly spaced along a circle.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The Path object to populate with coordinates
    spacing : float, optional
        Distance between sites along the path.
    N : int, optional
        Number of sites in the cyclic path.
    radius : float, optional
        The radius (nm) of the cyclic path.
    closed : bool, default True
        If `True` the cyclic path is closed by bonding the first and last sites together

    Notes
    -----
    Only two of spacing, N and radius can be defined, as the third
    is determined by the other two.
    """
    n_params = sum(1 for i in (spacing, N, radius) if i is not None)
    if n_params != 2:
        raise ValueError("You must specify only 2 of spacing, N and radius.")

    if spacing and N:
        radius = (N * spacing) / (2 * np.pi)
    elif radius and spacing:
        N = int((2 * np.pi * radius) / spacing)
    else:
        spacing = (2 * np.pi * radius) / N

    angles = np.arange(0, 2 * np.pi, (2 * np.pi) / N)
    path.coordinates = np.array(
        [(np.cos(a) * radius, np.sin(a) * radius, 0) for a in angles]
    )
    path.beads = [bead_name] * N
    path.form_linear_bond_graph()
    if closed:
        path.bond_graph.add_edge(0, N - 1)
        path.bond_graph.add_edge(0, N - 1)


def knot(path, spacing, N, m, closed=True, bead_name="_A"):
    """Generate a knot path.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The Path object to populate with coordinates
    spacing : float (nm)
        The spacing between sites along the path.
    N : int
        The number of total sites in the path.
    m : int in [3, 4, 5]
        The number of crossings in the knot.
        3 gives the trefoil knot, 4 gives the figure 8 knot and 5 gives the cinquefoil knot.
    closed : bool, default True
        If `True` the cyclic path is closed by bonding the first and last sites together
    """
    # Generate dense sites first, sample actual ones later from spacing
    t_dense = np.linspace(0, 2 * np.pi, 5000)

    # Trefoil knot
    if m == 3:
        R, r = 1.0, 0.3
        x = (R + r * np.cos(3 * t_dense)) * np.cos(2 * t_dense)
        y = (R + r * np.cos(3 * t_dense)) * np.sin(2 * t_dense)
        z = r * np.sin(3 * t_dense)
    # Figure-eight
    elif m == 4:
        x = (2 + np.cos(2 * t_dense)) * np.cos(3 * t_dense)
        y = (2 + np.cos(2 * t_dense)) * np.sin(3 * t_dense)
        z = np.sin(4 * t_dense)
    # Cinquefoil knot
    elif m == 5:
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
    L_target = (N - 1) * spacing

    # Scale to match target contour length
    scale = L_target / base_length
    coords_dense *= scale
    arc_lengths *= scale

    # Resample uniformly along arc length
    desired_arcs = np.linspace(0, L_target, N, endpoint=False)
    x_interp = interp1d(arc_lengths, coords_dense[:, 0])(desired_arcs)
    y_interp = interp1d(arc_lengths, coords_dense[:, 1])(desired_arcs)
    z_interp = interp1d(arc_lengths, coords_dense[:, 2])(desired_arcs)
    path.coordinates = np.stack((x_interp, y_interp, z_interp), axis=1)

    path.beads = [bead_name] * N
    path.form_linear_bond_graph()
    if closed:
        path.bond_graph.add_edge(0, len(path.coordinates) - 1)
        path.bond_graph.add_edge(0, len(path.coordinates) - 1)


def helix(
    path, N, radius, rise, twist, right_handed=True, bottom_up=True, bead_name="_A"
):
    """Generate helical path.

    Parameters:
    -----------
    path : mbuild.path.Path, required
        The Path object to populate with coordinates
    N : int, required
        Number of sites in the path
    radius : float, required
        Radius of the helix (nm)
    rise : float, required
        Rise per site on path (nm)
    twist : float, required
        Twist per site in path (degrees)
    right_handed : bool, default True
        Set the handedness of the helical twist
    bottom_up : bool, default True
        If True, the twist is in the positive Z direction
    """
    coordinates = np.zeros((N, 3))
    indices = reversed(range(N)) if not bottom_up else range(N)

    for i in indices:
        angle = np.deg2rad(i * twist)
        if not right_handed:
            angle *= -1
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        z = i * rise if bottom_up else -i * rise
        coordinates[i] = (x, y, z)

    path.coordinates = coordinates
    path.beads = [bead_name] * N
    path.form_linear_bond_graph()


def spiral_2D(path, N, a, b, spacing, bead_name="_A"):
    """Generate a 2D spiral path in the XY plane.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The Path object to populate with coordinates
    N : int, required
        Number of sites in the path
    a : float, required
        The initial radius (nm)
    b : float, required
        Determines how fast radius grows per angle increment (r = a + bθ)
    spacing : float, required
        Distance between adjacent sites (nm)
    """
    coordinates = np.zeros((N, 3))
    theta = 0.0

    for i in range(N):
        r = a + b * theta
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        z = 0.0
        coordinates[i] = (x, y, z)
        # Estimate next angle increment based on arc length
        ds_dtheta = np.sqrt((r) ** 2 + b**2)
        dtheta = spacing / ds_dtheta
        theta += dtheta

    path.coordinates = coordinates
    path.beads = [bead_name] * N
    path.form_linear_bond_graph()


def zigzag(
    path,
    N,
    spacing=1.0,
    angle_deg=120.0,
    sites_per_segment=4,
    plane="xy",
    bead_name="_A",
):
    """Generates a path following a zig-zag pattern in a given plane.

    Parameters
    ----------
    path : mbuild.path.Path, required
        The Path object to populate with coordinates
    N : int, required
        Number of sites in the path
    spacing : float, default = 1.0 nm
        The distance (nm) between consecutive sites along the path.
    angle_deg : float, default = 120.
        The rotation (degrees) applied between segments
    sites_per_segment : int, default = 4
        The number of sites before rotating and beginning next segment.
    plane : str, default = "xy"
        The plane that the sites in the path occupy
    """
    if N % sites_per_segment != 0:
        raise ValueError("N must be evenly divisible by sites_per_segment")

    angle_rad = np.deg2rad(angle_deg)
    direction = np.array([1.0, 0.0])
    position = np.zeros(2)
    coords_2d = []
    step_count = 0
    segment_count = 0

    for i in range(N):
        coords_2d.append(position.copy())
        position += spacing * direction
        step_count += 1
        # Rotate
        if step_count == sites_per_segment:
            step_count = 0
            sign = -1 if segment_count % 2 == 0 else 1
            rot_matrix = np.array(
                [
                    [np.cos(sign * angle_rad), -np.sin(sign * angle_rad)],
                    [np.sin(sign * angle_rad), np.cos(sign * angle_rad)],
                ]
            )
            rot_matrix = np.array(
                [
                    [np.cos(sign * angle_rad), -np.sin(sign * angle_rad)],
                    [np.sin(sign * angle_rad), np.cos(sign * angle_rad)],
                ]
            )
            direction = rot_matrix @ direction
            segment_count += 1

    # Map into 3D space based on chosen plane
    coordinates = np.zeros((N, 3))
    for i, (x2d, y2d) in enumerate(coords_2d):
        if plane == "xy":
            coordinates[i] = (x2d, y2d, 0)
        elif plane == "xz":
            coordinates[i] = (x2d, 0, y2d)
        elif plane == "yz":
            coordinates[i] = (0, x2d, y2d)

    path.coordinates = coordinates
    path.beads = [bead_name] * N  # TODO: Should extend this, not assign it
    path.beads = [bead_name] * N  # TODO: Should extend this, not assign it
    path.form_linear_bond_graph()


def hard_sphere_random_walk(
    path=None,
    bead_name="_A",
    bond_length=0.15,
    radius=0.1,
    rw_angles=None,
    termination=None,
    volume_constraint=None,
    bias=None,
    connectivity="linear",
    initial_point=None,
    seed=42,
    trial_batch_size=20,
    tolerance=1e-5,
    chunk_size=512,
    run_on_gpu=False,
):
    """Generates coordinates from a self avoiding random walk using
    fixed bond lengths, hard spheres, and minimum and maximum angles
    formed by 3 consecutive points.

    Parameters:
    -----------
    path : mbuild.path.Path, default None.
        The Path object to populate with coordinates. Creates a new path object if not passed.
    bead_name : str, default `'_A'`
        Name to assign to nodes added during this walk.
    bond_length : float, default 0.15 nm
        Fixed bond length between 2 coordinates.
    radius : float, default 0.1 nm
        Radius of sites used in checking for overlaps.
    rw_angles : tuple or dict or np.array or AnglesSampler, default None
        Set the angle sampling method. A tuple of (min_val, max_val) sets the uniform distribution.
        The default value of None sets the uniform angle sampling from (np.pi/2, np.pi). Can also
        use a Gaussian distribution by passing a dict with keys {'loc':mean, 'scale':std}.
        Finally, a numpy array of 1D or 2D array of numpy values can be passed, which will be sampled
        via numpy.random.choice method. The 2D case provides a set of weights.
    termination : termination condition, required
        Termination condition for the random walk. If an integer is passed,
        will terminate after reaching that number of sites. Can also pass a tuple of
        will terminate after reaching that number of sites. Can also pass a tuple of
        termination conditions, or a mbuild.path.termination.Termination object.
    volume_constraint : mbuild.path.constraints.Constraint, optional
        Used to reject moves which are outside of the volume constraint.
        See mbuild.utils.Constraint objects.
    bias : mbuild.path.bias.Bias class,
        Causes the selected points to be biased towards some criterion.
    connectivity : str, default "linear"
        Will be used to connect the bond_graph of the walk post generation based on a specified method.
        See path._connect_edges for different options.
    initial_point : array-like or int, optional
        Used as the coordinate for the first site in this random walk path. If an integer is
        Used as the coordinate for the first site in this random walk path. If an integer is
        passed, look in coordinates of passed path object, and grab the starting coordiantes from there.
    seed : int, default = 42
        Random seed
    trial_batch_size : int, default = 20
        The number of trial moves to attempt in parallel for each step.
    tolerance : float, default = 1e-5
        Tolerance used for rounding and checking for overlaps.
    chunk_size : int, default = 512
        Size of coordinate chunks to allocate
    run_on_gpu : bool, default = False
        If True and CUDA path utilities are available, use GPU-accelerated
        implementations.
    """
    # Create state object to track random walk progress
    state = RandomWalkState(
        bond_length=bond_length,
        radius=radius,
        angles_sampler=rw_angles,
        bead_name=bead_name,
        initial_point=initial_point,
        previous_count=len(path.coordinates) if path else 0,
        connectivity=connectivity,
        seed=seed,
        volume_constraint=volume_constraint,
        tolerance=tolerance,
        trial_batch_size=int(trial_batch_size),
        chunk_size=chunk_size,
        run_on_gpu=bool(run_on_gpu) and _get_cuda_available(),
    )

    if path is None:  # Create empty path
        path = Path()

    if termination is None:  # handle default termination args
        raise ValueError(
            "Please pass viable termination to hard_sphere_random_walk from `mbuild.path.termination.Termination`"
        )
    elif isinstance(termination, int):
        termination = Termination(NumSites(termination))
    elif isinstance(termination, (tuple, list)):
        termination = Termination(termination)
    elif isinstance(termination, Terminator):
        termination = Termination(termination)
    elif not isinstance(termination, Termination):
        raise ValueError(f"Bad input {termination=} to hard_sphere_random_walk")
    state.termination = termination
    state.termination._attach_path(path, state)

    # Create RNG state
    rng = np.random.default_rng(seed + len(path.coordinates))
    state.rng = rng

    # Set up PBC info from volume constraints
    if isinstance(volume_constraint, CuboidConstraint):
        pbc = volume_constraint.pbc
        box_lengths = volume_constraint.box_lengths.astype(np.float32)
    elif isinstance(volume_constraint, CylinderConstraint):
        pbc = (False, False, volume_constraint.periodic_height)
        box_lengths = np.array(
            [
                volume_constraint.radius * 2,
                volume_constraint.radius * 2,
                volume_constraint.height,
            ]
        ).astype(np.float32)
    else:
        pbc = (None, None, None)
        box_lengths = (None, None, None)

    # Set up bias conditions
    if bias:
        bias._attach_path(path, state)
        state.bias = bias

    # Initialize coordinates based on starting conditions
    if not path.coordinates.size == 0:
        coordinates = np.concatenate(
            (
                path.coordinates.astype(np.float32),
                np.zeros((chunk_size, 3), dtype=np.float32),
            ),
            axis=0,
        )
        state.count = len(path.coordinates)  # starting index
        state.count = len(path.coordinates)  # starting index
    else:
        coordinates = np.zeros((chunk_size, 3), dtype=np.float32)
        state.count = 0

    state.init_count = state.count

    # Set start time for wall time terminator
    state.start_time = time.time()

    # Select methods for random walk
    if state.run_on_gpu:
        from mbuild.path.path_utils_gpu import check_path_split

        logger.info("Running hard_sphere_random_walk on a CUDA device.")
        check_path_gpu = check_path_split
    else:
        check_path_gpu = None

    check_path_cpu = check_path
    next_step = random_coordinate

    # Set the first two coordinates
    for num_tries in range(100):
        initial_xyz = get_initial_point(
            state, coordinates[: state.count], check_path_cpu, next_step
        )

        if state.check_termination(path, coordinates):
            return path

        # Get the second coordinate
        second_xyz = get_second_point(
            state, coordinates[: state.count], check_path_cpu, initial_xyz
        )
        if second_xyz is not None:
            break
        elif num_tries == 99:
            raise PathConvergenceError(
                f"Failed after {num_tries + 1} to generate a starting point. System is probably too densely packed."
            )
        elif state.initial_point is not None:
            raise PathConvergenceError(
                f"Failed to initiate random walk with {initial_point=}. Try a different initial_point."
            )

    # print(f"{initial_xyz=}, {second_xyz=}")
    coordinates[state.count] = initial_xyz
    state.count += 1
    coordinates[state.count] = second_xyz

    if state.check_termination(path, coordinates):
        return path

    # Prepare GPU static points if using GPU
    if state.run_on_gpu:
        from numba import cuda

        static_parts = []
        if state.init_count > 0:
            static_parts.append(coordinates[: state.init_count])
        if static_parts:
            static_points = np.concatenate(static_parts).astype(np.float32)
            state.gpu_static_points = cuda.to_device(static_points)

    # Main random walk loop
    walk_finished = False
    while not walk_finished:
        batch_angles, batch_vectors = generate_trials(state)
        candidates = next_step(
            pos1=coordinates[state.count],
            pos2=coordinates[state.count - 1],
            bond_length=bond_length,
            thetas=batch_angles,
            r_vectors=batch_vectors,
        )

        if state.volume_constraint:  # should allow for pbc
            is_inside_mask = volume_constraint.is_inside(
                points=candidates, buffer=radius
            )
            candidates = candidates[is_inside_mask]

        if state.bias:
            candidates = bias(candidates=candidates)

        existing_points = coordinates[: state.count + 1]

        if any(pbc):
            candidates = volume_constraint.mins + np.mod(
                candidates - volume_constraint.mins, box_lengths
            )

        accept_xyz = None
        if state.run_on_gpu and len(candidates) > 0:
            dynamic_points = coordinates[state.init_count : state.count + 1]
            valid_mask = check_path_gpu(
                state.gpu_static_points,
                dynamic_points,
                candidates,
                radius,
                tolerance,
            )
            valid_candidates = candidates[valid_mask]
            if len(valid_candidates) > 0:
                accept_xyz = valid_candidates[0]
        else:
            for xyz in candidates:
                if check_path_cpu(
                    existing_points=existing_points,
                    new_point=xyz,
                    radius=radius,
                    tolerance=tolerance,
                ):
                    accept_xyz = xyz
                    break

        if accept_xyz is not None:
            coordinates[state.count + 1] = accept_xyz
            state.count += 1

        state.attempts += 1

        # Extend coordinates array if we're running out of space
        if state.count + 2 >= len(coordinates):
            path.coordinates = coordinates  # Save progress first
            path._extend_coordinates(N=chunk_size)
            coordinates = path.coordinates

        walk_finished = termination.is_met()

    state.check_termination(path, coordinates)

    return path


class RandomWalkState:
    """Tracks state and configuration for a hard_sphere_random_walk.


    This class encapsulates all the bookkeeping information needed during
    a random walk, keeping the Path object clean of implementation details.


    Attributes
    ----------
    bond_length : float
        Fixed bond length between coordinates
    radius : float
        Radius of sites for overlap checking
    min_angle : float
        Minimum angle (radians) for selecting next step
    max_angle : float
        Maximum angle (radians) for selecting next step
    count : int
        Current number of successfully placed sites
    init_count : int
        Number of sites at the start of this random walk (from previous path)
    attempts : int
        Total number of attempted moves
    start_time : float
        Time when the random walk started (for WallTime terminator)
    initial_point : np.ndarray or None
        Specified initial coordinate
    seed : int
        Random seed
    volume_constraint : Constraint or None
        Volume constraint for the walk
    tolerance : float
        Tolerance for overlap checking
    trial_batch_size : int
        Number of trial moves per step
    chunk_size : int
        Size of coordinate chunks to allocate
    run_on_gpu : bool
        Whether GPU acceleration is being used
    gpu_static_points : device array or None
        GPU array of static points for overlap checking
    """

    def __init__(
        self,
        bond_length,
        radius,
        angles_sampler,
        bead_name,
        initial_point=None,
        previous_count=0,
        connectivity=None,
        seed=42,
        volume_constraint=None,
        termination=None,
        bias=None,
        tolerance=1e-5,
        trial_batch_size=20,
        chunk_size=512,
        run_on_gpu=False,
    ):
        self.bond_length = bond_length
        self.radius = radius
        if bond_length < radius:
            raise ValueError(
                "Bond length should be greater than radius to prevent overlaps."
            )
        if angles_sampler is None:
            self.angles = AnglesSampler(
                "uniform", {"low": np.pi / 2, "high": np.pi}, seed
            )
        elif isinstance(angles_sampler, tuple):
            self.angles = AnglesSampler(
                "uniform", {"low": angles_sampler[0], "high": angles_sampler[1]}, seed
            )
        elif (
            isinstance(angles_sampler, dict)
            and angles_sampler.get("loc")
            and angles_sampler.get("scale")
        ):
            self.angles = AnglesSampler("normal", angles_sampler, seed)
        elif isinstance(angles_sampler, np.ndarray):
            if angles_sampler.ndim == 1:
                kwargs = {"a": angles_sampler}
            elif angles_sampler.ndim == 2:
                kwargs = {"a": angles_sampler[0], "p": angles_sampler[1]}
            self.angles = AnglesSampler("choice", kwargs, seed)
        elif isinstance(angles_sampler, AnglesSampler):
            self.angles = angles_sampler
        else:
            raise ValueError(
                f"Please provide a reasonable value to set the rw_angles. Passed {angles_sampler}"
            )
        self.bead_name = bead_name
        if hasattr(initial_point, "__len__") and len(initial_point) == 3:
            self.initial_point = np.asarray(initial_point)
        else:
            self.initial_point = initial_point
        self.previous_count = previous_count
        self.connectivity = connectivity
        self.seed = seed
        self.volume_constraint = volume_constraint
        self.termination = termination
        self.tolerance = tolerance
        self.bias = bias
        self.trial_batch_size = trial_batch_size
        self.chunk_size = chunk_size
        self.run_on_gpu = run_on_gpu

        # State tracking
        self.count = 0
        self.init_count = 0
        self.attempts = 0
        self.start_time = None
        self.gpu_static_points = None

    def check_termination(self, path, coordinates):
        """Examine and process termination if we have reached."""
        if self.termination.is_met():
            if self.termination.success:
                logger.info("Random walk successful.")
                path.coordinates = coordinates[: self.count + 1]
            else:
                logger.warning("Random walk not successful.")
                logger.warning(self.termination.summarize())
            self.termination._clean()
            if self.bias:
                self.bias._clean()
            path._extend_bond_graph(self.bead_name)
            if isinstance(
                self.initial_point, int
            ):  # make sure to build from previous point instead of last point
                path._connect_edges(
                    self.connectivity,
                    np.arange(self.previous_count, self.count + 1),
                    self.initial_point,
                )
            else:  # build bond graph, and connect to last index in previous path coordinates
                path._connect_edges(
                    self.connectivity,
                    np.arange(self.previous_count, self.count + 1),
                    self.previous_count,
                )
            if isinstance(
                self.initial_point, int
            ):  # make sure to build from previous point instead of last point
                path._connect_edges(
                    self.connectivity,
                    np.arange(self.previous_count, self.count + 1),
                    self.initial_point,
                )
            else:  # build bond graph, and connect to last index in previous path coordinates
                path._connect_edges(
                    self.connectivity,
                    np.arange(self.previous_count, self.count + 1),
                    self.previous_count,
                )
            path._extend_beads(self.bead_name)
            return True
        return False


def crosslink(
    path,
    bead_name="_R",
    backbone_name="_A",
    radius=0.1,
    excluded_bond_depth=2,
    n_connection_sites=2,
    volume_constraint=None,
    initial_point=None,
    seed=42,
    chunk_size=512,
    run_on_gpu=False,
):
    """
    Create a crosslink node that bonds to n_connection_sites backbone beads.

    Adds a new node with bead_name to path.bond_graph, positioned near
    and bonded to n_connection_sites backbone beads within the specified radius.

    Parameters
    ----------
    path : Path
        The Path object containing coordinates and bond_graph
    bead_name : str, default "_R"
        Name for the crosslink bead
    backbone_name : str, default "_A"
        Name for backbone beads to search for
    radius : float, default 0.1
        Search radius for finding nearby backbone beads
    n_connection_sites : int, default 2
        Number of backbone beads to bond to
    initial_point : int or array-like, optional
        Starting point (node index or xyz coordinate) to search around
    seed : int, default 42
        Random seed for reproducibility
    chunk_size : int, default 512
        Chunk size for batch processing (used if extending coordinates)
    run_on_gpu : bool, default False
        Whether to use GPU acceleration via numba

    Returns
    -------
    Path
        The modified path object with the new crosslink node
    """
    rng = np.random.default_rng(seed + len(path.coordinates))

    # Find all backbone beads
    backbone_nodes = [
        node
        for node in path.bond_graph.nodes()
        if (
            path.beads[node] == backbone_name
            # and path.bond_graph.degree[node] <= 2
        )  # TODO: Multiple crosslink sites on one backbone too
    ]  # value references global node index
    backbone_subgraph = path.bond_graph.subgraph(
        backbone_nodes
    )  # TODO: make a path function?

    if len(backbone_nodes) == 0:
        raise ValueError(f"No backbone beads with name '{backbone_name}' found in path")

    if len(backbone_nodes) < n_connection_sites:
        raise ValueError(
            f"Not enough backbone beads ({len(backbone_nodes)}) for "
            f"{n_connection_sites} connection sites"
        )

        # Set up PBC info from volume constraints
    if isinstance(volume_constraint, CuboidConstraint):
        pbc = volume_constraint.pbc
        box_lengths = volume_constraint.box_lengths.astype(np.float32)
    elif isinstance(volume_constraint, CylinderConstraint):
        pbc = (False, False, volume_constraint.periodic_height)
        box_lengths = np.array(
            [
                volume_constraint.radius * 2,
                volume_constraint.radius * 2,
                volume_constraint.height,
            ]
        ).astype(np.float32)
    else:
        pbc = np.array([False, False, False], dtype=bool)
        box_lengths = np.array([np.inf, np.inf, np.inf], dtype=np.float32)

    # Get coordinates of all backbone nodes
    candidate_nodes = [
        node for node in backbone_nodes if path.bond_graph.degree[node] <= 2
    ]
    candidate_coords = np.array(path.coordinates[candidate_nodes], dtype=np.float32)

    # get reference points
    def get_reference_points(path, initial_point):
        """Create reference points for finding candidates.

        Returns
        -------
        nodesArray: np.array
            index of global bond_graph nodes that are viable starting points -> [0,2,10...]
        coordsArray: np.array
            each value matches path.coordinates[nodesList]
        """
        if initial_point is not None:
            if isinstance(initial_point, (int, np.integer)):
                # Use coordinate of specified node
                if initial_point not in path.bond_graph.nodes:
                    raise ValueError(f"Node {initial_point} not found in bond_graph")
                nodesArray = np.array([initial_point])
                coordsArray = np.array([path.coordinates[initial_point]])
            else:
                # Use provided coordinate
                initial_point32 = np.asarray(initial_point, dtype=np.float32)
                sq_distances = calculate_sq_distances(
                    initial_point32, candidate_coords, pbc=pbc, box_lengths=box_lengths
                )
                nodesArray = np.argsort(sq_distances)
                coordsArray = path.coordinates[nodesArray]
        else:
            # Randomly select a backbone node as reference
            nodesArray = rng.choice(
                candidate_nodes, size=len(candidate_nodes), replace=False
            )
            coordsArray = path.coordinates[nodesArray]

        return nodesArray, coordsArray

    ref_nodes, ref_coords = get_reference_points(path, initial_point)

    found_ref = False  # flag to check all ref_nodes
    for ref_node, ref_coord in zip(ref_nodes, ref_coords):
        selected_nodes = [ref_node]  # first choice is ref
        # GPU-accelerated distance calculation
        sq_distances = calculate_sq_distances(
            ref_coord, candidate_coords, pbc=pbc, box_lengths=box_lengths
        )
        distances = np.sqrt(sq_distances)

        # Find candidates within radius
        within_radius_mask = distances <= (radius) * 2  # twice radiu
        possible_pairs = np.where(within_radius_mask)[0]
        # Verify starting point or return early
        if len(possible_pairs) < n_connection_sites - 1:
            continue

        closest_paired_nodes = possible_pairs[np.argsort(distances[possible_pairs])]
        excluded_nodes = set(
            nx.single_source_shortest_path_length(
                backbone_subgraph, ref_node, cutoff=excluded_bond_depth
            ).keys()
        )
        # import pdb; pdb.set_trace()
        for idx in closest_paired_nodes:
            node = candidate_nodes[idx]  # temp replace
            # node = possible_pairs[idx] # is index == value ??
            if node in excluded_nodes:
                continue

            selected_nodes.append(int(node))

            # Stop if we have enough connection sites
            if len(selected_nodes) >= n_connection_sites:
                found_ref = True
                break  # break twice
        if found_ref:
            break

    # Verify enough final viable crosslink
    if not found_ref:
        n_clinks = sum([bead == bead_name for bead in path.beads])
        raise PathConvergenceError(
            f"Only found {len(selected_nodes)} non-neighboring backbone beads "
            f"within radius {radius}, need {n_connection_sites}."
            f"\nMaximum crossinks found are {n_clinks}. "
            "Ways to increase crosslinking:\nIncrease radius"
            "\nPack at higher density\nRelax structure."
        )

    # Calculate position for new crosslink node (centroid of selected beads)
    selected_coords = np.array(path.coordinates[selected_nodes])
    crosslink_position = np.mean(selected_coords, axis=0)

    # Add new node to path
    path.append_coordinates(crosslink_position, bead_name)
    new_node_idx = len(path.coordinates) - 1  # add as last index

    # Add edges from crosslink node to selected backbone nodes
    for backbone_node in selected_nodes:
        path.bond_graph.add_edge(
            int(new_node_idx),
            int(backbone_node),
            bond_type=(bead_name, backbone_name),
        )

    return path


class CrosslinkWalkState:
    # TODO
    pass


_CUDA_AVAILABLE = None


def _get_cuda_available():
    """Check if numba can access CUDA runtime."""
    global _CUDA_AVAILABLE
    if _CUDA_AVAILABLE is None:
        try:
            from numba import cuda

            _CUDA_AVAILABLE = cuda.is_available()
        except Exception:
            _CUDA_AVAILABLE = False
    return _CUDA_AVAILABLE
