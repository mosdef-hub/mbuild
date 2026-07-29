"""Crosslinking module for polymer path structures.

Handles placement of crosslinker geometries between backbone beads,
with proper overlap avoidance and bond length enforcement under PBC.
"""

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

from mbuild.exceptions import PathConvergenceError
from mbuild.path.build import Path
from mbuild.path.constraints import CuboidConstraint, CylinderConstraint

# =============================================================================
# CrosslinkerGeometry
# =============================================================================


class CrosslinkerGeometry(Path):
    """Rigid crosslinker geometry with connection site metadata.

    Parameters
    ----------
    coordinates : array-like, shape (N, 3)
        Positions of beads in the crosslinker geometry.
    bond_graph : networkx.Graph, optional
        Internal bond topology. If None, creates a graph with no edges.
    bead_name : str or array-like of str
        Name(s) for each bead. Strings are prefixed with '_' if not already.
    connection_sites : list of int
        Node indices that bond to external backbone beads. May contain
        duplicates to indicate multiple connections from the same physical node.
    """

    def __init__(
        self, coordinates, bond_graph=None, bead_name="_R", connection_sites=None
    ):
        coordinates = np.asarray(coordinates, dtype=np.float32)
        self.coordinates = coordinates - coordinates.mean(axis=0)

        if bond_graph is None:
            bond_graph = nx.Graph()
            bond_graph.add_nodes_from(range(self.n_sites))
        self.bond_graph = bond_graph

        if isinstance(bead_name, str):
            if not bead_name.startswith("_"):
                bead_name = "_" + bead_name
            self.beads = np.full(self.n_sites, bead_name, dtype="U10")
        else:
            self.beads = np.asarray(bead_name, dtype="U10")

        if connection_sites is None:
            connection_sites = list(range(self.n_sites))
        self.connection_sites = list(connection_sites)

        self.validate_connection_sites()

    @property
    def n_sites(self):
        """Number of beads in the crosslinker."""
        return len(self.coordinates)

    @property
    def n_connections(self):
        """Total number of connections (including duplicates)."""
        return len(self.connection_sites)

    @property
    def unique_connection_sites(self):
        """Unique physical node indices used as connection sites."""
        return list(dict.fromkeys(self.connection_sites))

    @property
    def internal_bonds(self):
        """List of (u, v) edges within the crosslinker."""
        return list(self.bond_graph.edges())

    def validate_connection_sites(self):
        """Raise if any connection site index exceeds the number of sites."""
        largest_index = max(self.connection_sites)
        if largest_index >= self.n_sites:
            raise PathConvergenceError(
                f"connection site {largest_index} is too large for {self.n_sites}"
            )

    def copy(self):
        """Return a deep copy of this crosslinker geometry."""
        return CrosslinkerGeometry(
            coordinates=self.coordinates.copy(),
            bond_graph=self.bond_graph.copy(),
            bead_name=self.beads.copy(),
            connection_sites=list(self.connection_sites),
        )

    def recenter(self):
        """Shift coordinates so the centroid is at the origin."""
        if len(self.coordinates) > 0:
            self.coordinates -= self.coordinates.mean(axis=0)

    @classmethod
    def from_path(cls, path, connection_sites):
        """Create a CrosslinkerGeometry from an existing Path object.

        Parameters
        ----------
        path : Path
            Source path providing coordinates, bond_graph, and bead names.
        connection_sites : list of int
            Which nodes serve as external connection points.
        """
        return cls(
            coordinates=path.coordinates.copy(),
            bond_graph=path.bond_graph.copy(),
            bead_name=path.beads.copy(),
            connection_sites=connection_sites,
        )

    @classmethod
    def single_site(cls, bead_name="_R", n_connections=2):
        """A single bead that bonds to n_connections backbone beads.

        Parameters
        ----------
        bead_name : str
            Name for the crosslinker bead.
        n_connections : int
            Number of backbone beads this single site connects to.
        """
        return cls(
            coordinates=np.zeros((1, 3), dtype=np.float32),
            bead_name=bead_name,
            connection_sites=[0] * n_connections,
        )

    @classmethod
    def equilateral_triangle(
        cls, bond_length=0.27, bead_name="_R", connection_sites=None
    ):
        """Three beads arranged in an equilateral triangle.

        Parameters
        ----------
        bond_length : float
            Edge length of the triangle.
        bead_name : str
            Name for each bead.
        connection_sites : list of int, optional
            Defaults to [0, 1, 2].
        """
        R = bond_length / np.sqrt(3)
        angles = np.array([0, 2 * np.pi / 3, 4 * np.pi / 3])
        positions = np.column_stack(
            [R * np.cos(angles), R * np.sin(angles), np.zeros(3)]
        ).astype(np.float32)
        G = nx.cycle_graph(3)
        if connection_sites is None:
            connection_sites = [0, 1, 2]
        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_name,
            connection_sites=connection_sites,
        )

    @classmethod
    def linear(cls, n_sites=2, bond_length=0.27, bead_name="_R", connection_sites=None):
        """A linear chain of beads along the x-axis.

        Parameters
        ----------
        n_sites : int
            Number of beads in the chain.
        bond_length : float
            Spacing between adjacent beads.
        bead_name : str
            Name for each bead.
        connection_sites : list of int, optional
            Defaults to the two endpoints [0, n_sites - 1].
        """
        positions = np.zeros((n_sites, 3), dtype=np.float32)
        positions[:, 0] = np.arange(n_sites) * bond_length
        G = nx.path_graph(n_sites)
        if connection_sites is None:
            connection_sites = [0, n_sites - 1]
        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_name,
            connection_sites=connection_sites,
        )

    @classmethod
    def square(cls, bond_length=0.27, bead_name="_R", connection_sites=None):
        """Four beads at the corners of a square.

        Parameters
        ----------
        bond_length : float
            Edge length of the square.
        bead_name : str
            Name for each bead.
        connection_sites : list of int, optional
            Defaults to [0, 1, 2, 3].
        """
        half = bond_length / 2.0
        positions = np.array(
            [[half, half, 0], [-half, half, 0], [-half, -half, 0], [half, -half, 0]],
            dtype=np.float32,
        )
        G = nx.cycle_graph(4)
        if connection_sites is None:
            connection_sites = [0, 1, 2, 3]
        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_name,
            connection_sites=connection_sites,
        )

    @classmethod
    def tetrahedral(cls, bond_length=0.27, bead_name="_R", connection_sites=None):
        """Four beads at the vertices of a regular tetrahedron.

        Parameters
        ----------
        bond_length : float
            Edge length of the tetrahedron.
        bead_name : str
            Name for each bead.
        connection_sites : list of int, optional
            Defaults to [0, 1, 2, 3].
        """
        positions = np.array(
            [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=np.float32
        )
        positions *= bond_length / np.linalg.norm(positions[0] - positions[1])
        G = nx.complete_graph(4)
        if connection_sites is None:
            connection_sites = [0, 1, 2, 3]
        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_name,
            connection_sites=connection_sites,
        )

    @classmethod
    def trigonal_bipyramidal(
        cls, bond_length=0.27, bead_name="_R", connection_sites=None
    ):
        """Five beads in a trigonal bipyramidal arrangement.

        Three equatorial beads plus two axial beads. All equatorial-equatorial
        and equatorial-axial edges are present.

        Parameters
        ----------
        bond_length : float
            Equatorial edge length.
        bead_name : str
            Name for each bead.
        connection_sites : list of int, optional
            Defaults to [0, 1, 2, 3, 4].
        """
        R_eq = bond_length / np.sqrt(3)
        eq_angles = np.array([0, 2 * np.pi / 3, 4 * np.pi / 3])
        equatorial = np.column_stack(
            [R_eq * np.cos(eq_angles), R_eq * np.sin(eq_angles), np.zeros(3)]
        )
        axial_dist = bond_length * np.sqrt(2.0 / 3.0)
        axial = np.array([[0, 0, axial_dist], [0, 0, -axial_dist]])
        positions = np.vstack([equatorial, axial]).astype(np.float32)

        G = nx.Graph()
        G.add_nodes_from(range(5))
        G.add_edges_from([(0, 1), (1, 2), (0, 2)])
        for eq in range(3):
            G.add_edges_from([(eq, 3), (eq, 4)])
        if connection_sites is None:
            connection_sites = [0, 1, 2, 3, 4]
        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_name,
            connection_sites=connection_sites,
        )

    @classmethod
    def pentagon(cls, bond_length=0.27, bead_name="_R", connection_sites=None):
        """Five beads arranged in a regular pentagon.

        Parameters
        ----------
        bond_length : float
            Edge length of the pentagon.
        bead_name : str
            Name for each bead.
        connection_sites : list of int, optional
            Defaults to [0, 1, 2, 3, 4].
        """
        R = bond_length / (2 * np.sin(np.pi / 5))
        angles = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        positions = np.column_stack(
            [R * np.cos(angles), R * np.sin(angles), np.zeros(5)]
        ).astype(np.float32)
        G = nx.cycle_graph(5)
        if connection_sites is None:
            connection_sites = [0, 1, 2, 3, 4]
        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_name,
            connection_sites=connection_sites,
        )

    @classmethod
    def from_edges(cls, coordinates, edges, bead_name="_R", connection_sites=None):
        """Create a crosslinker from explicit coordinates and edge list.

        Parameters
        ----------
        coordinates : array-like, shape (N, 3)
            Bead positions.
        edges : list of (int, int)
            Pairs of node indices defining internal bonds.
        bead_name : str
            Name for each bead.
        connection_sites : list of int, optional
            Defaults to all nodes.
        """
        coordinates = np.asarray(coordinates, dtype=np.float32)
        n = len(coordinates)
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(edges)
        if connection_sites is None:
            connection_sites = list(range(n))
        return cls(
            coordinates=coordinates,
            bond_graph=G,
            bead_name=bead_name,
            connection_sites=connection_sites,
        )


# =============================================================================
# PBC Utilities
# =============================================================================


def _pbc_delta(a, b, box_lengths, pbc):
    """Minimum-image displacement vector (a - b), vectorized.

    Parameters
    ----------
    a, b : array-like, shape (..., 3)
        Position(s).
    box_lengths : array-like, shape (3,)
        Box dimensions.
    pbc : array-like of bool, shape (3,)
        Which dimensions are periodic.

    Returns
    -------
    delta : ndarray, shape (..., 3)
    """
    delta = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    for dim in range(3):
        if pbc[dim] and box_lengths[dim] > 0:
            delta[..., dim] -= (
                np.rint(delta[..., dim] / box_lengths[dim]) * box_lengths[dim]
            )
    return delta


def _pbc_distance(a, b, box_lengths, pbc):
    """Minimum-image scalar distance between two points."""
    delta = _pbc_delta(a, b, box_lengths, pbc)
    return float(np.linalg.norm(delta, axis=-1))


def _get_pbc_info(volume_constraint):
    """Extract PBC flags and box lengths from a volume constraint.

    Returns
    -------
    pbc : ndarray of bool, shape (3,)
    box_lengths : ndarray of float, shape (3,)
    """
    if volume_constraint is None:
        return np.zeros(3, dtype=bool), np.zeros(3, dtype=np.float64)

    if isinstance(volume_constraint, CuboidConstraint):
        return (
            np.asarray(volume_constraint.pbc, dtype=bool),
            np.asarray(volume_constraint.box_lengths, dtype=np.float64),
        )
    elif isinstance(volume_constraint, CylinderConstraint):
        pbc = np.array(
            [False, False, getattr(volume_constraint, "pbc_z", False)], dtype=bool
        )
        box_lengths = np.array(
            [0.0, 0.0, getattr(volume_constraint, "length", 0.0)], dtype=np.float64
        )
        return pbc, box_lengths
    return np.zeros(3, dtype=bool), np.zeros(3, dtype=np.float64)


# =============================================================================
# Backbone Specification Parsing
# =============================================================================


def _parse_backbone_spec(backbone_name, connection_sites):
    """Parse backbone_name into per-connection-site specifications.

    Parameters
    ----------
    backbone_name : str or list
        Filter specification for eligible backbone beads.
    connection_sites : list
        The crosslinker's connection_sites list.

    Returns
    -------
    list : One spec entry per connection site. Each entry is a string
        (single bead filter) or tuple of strings (multi-bead filter).
    """
    n = len(connection_sites)
    if backbone_name is None:
        return None
    if isinstance(backbone_name, str):
        return [backbone_name] * n
    if len(backbone_name) != n:
        raise PathConvergenceError(
            f"backbone_name list length {len(backbone_name)} != "
            f"number of connection sites {n}"
        )
    return list(backbone_name)


def _beads_per_site(backbone_specs):
    """Return the number of beads each connection site bonds to.

    Parameters
    ----------
    backbone_specs : list
        Per-site specifications from _parse_backbone_spec.

    Returns
    -------
    list of int
    """
    if backbone_specs is None:
        return None
    counts = []
    for spec in backbone_specs:
        if (
            isinstance(spec, (list, tuple))
            and len(spec) > 0
            and isinstance(spec[0], str)
        ):
            counts.append(len(spec))
        else:
            counts.append(1)
    return counts


def _bead_matches_spec(bead_name, spec):
    """Check if a bead name matches a filter specification.

    Matching is underscore-prefix-insensitive: "_B" matches "B" and vice versa.

    Parameters
    ----------
    bead_name : str
        The bead's name.
    spec : str, list/tuple of str, or None
        None matches anything. A tuple matches if any element matches.

    Returns
    -------
    bool
    """
    if spec is None:
        return True
    if isinstance(spec, str):
        return (
            bead_name == spec
            or bead_name == "_" + spec
            or bead_name.lstrip("_") == spec.lstrip("_")
        )
    if isinstance(spec, (list, tuple)):
        return any(_bead_matches_spec(bead_name, s) for s in spec)
    return False


# =============================================================================
# Geometric Solvers
# =============================================================================


def _find_equidistant_point(
    bead_positions, target_distance, pbc=None, box_lengths=None
):
    """Find a point at target_distance from all given bead positions.

    Uses analytical solutions for n <= 3 beads and iterative projection
    for larger sets.

    Parameters
    ----------
    bead_positions : array-like, shape (N, 3)
    target_distance : float
    pbc : array-like of bool, optional (unused, reserved for future)
    box_lengths : array-like, optional (unused, reserved for future)

    Returns
    -------
    point : ndarray, shape (3,)
    """
    bead_positions = np.asarray(bead_positions, dtype=np.float64)
    n = len(bead_positions)

    if n == 0:
        return np.zeros(3, dtype=np.float64)
    if n == 1:
        return bead_positions[0] + np.array([0, 0, target_distance], dtype=np.float64)

    if n == 2:
        midpoint = (bead_positions[0] + bead_positions[1]) / 2.0
        half_sep = np.linalg.norm(bead_positions[1] - bead_positions[0]) / 2.0
        r_sq = target_distance**2 - half_sep**2
        if r_sq <= 0:
            return midpoint
        axis = bead_positions[1] - bead_positions[0]
        axis_norm = axis / (np.linalg.norm(axis) + 1e-30)
        perp = np.cross(axis_norm, np.array([0, 0, 1], dtype=np.float64))
        if np.linalg.norm(perp) < 1e-10:
            perp = np.cross(axis_norm, np.array([0, 1, 0], dtype=np.float64))
        perp /= np.linalg.norm(perp) + 1e-30
        return midpoint + perp * np.sqrt(r_sq)

    # General case (n >= 3)
    centroid = bead_positions.mean(axis=0)
    distances_from_centroid = np.linalg.norm(bead_positions - centroid, axis=1)

    if np.allclose(distances_from_centroid, target_distance, rtol=1e-6):
        return centroid

    # Symmetric arrangement: solution on axis through centroid perpendicular to bead plane
    if (
        np.std(distances_from_centroid) / (np.mean(distances_from_centroid) + 1e-30)
        < 1e-6
    ):
        r_centroid = distances_from_centroid[0]
        h_sq = target_distance**2 - r_centroid**2
        if h_sq <= 0:
            return centroid
        centered = bead_positions - centroid
        _, _, Vh = np.linalg.svd(centered)
        normal = Vh[-1]
        return centroid + normal * np.sqrt(h_sq)

    # Non-symmetric: iterative projection
    point = centroid.copy()
    for _ in range(100):
        shift = np.zeros(3, dtype=np.float64)
        for i in range(n):
            delta = point - bead_positions[i]
            dist = np.linalg.norm(delta)
            if dist < 1e-12:
                delta = np.array([1e-6, 1e-6, 1e-6], dtype=np.float64)
                dist = np.linalg.norm(delta)
            desired = bead_positions[i] + delta * (target_distance / dist)
            shift += desired - point
        shift /= n
        point += shift
        if np.linalg.norm(shift) < target_distance * 1e-10:
            break

    return point


def _can_reach_all_beads(
    bead_positions, crosslink_bond_length, tolerance, pbc, box_lengths
):
    """Check if a single point can be at crosslink_bond_length from all beads.

    Parameters
    ----------
    bead_positions : array-like, shape (N, 3)
    crosslink_bond_length : float
    tolerance : float
        Fractional tolerance on distance.
    pbc : array-like of bool, shape (3,)
    box_lengths : array-like, shape (3,)

    Returns
    -------
    feasible : bool
    target_point : ndarray or None
    """
    bead_positions = np.asarray(bead_positions, dtype=np.float64)
    n = len(bead_positions)
    r = crosslink_bond_length

    if n == 0:
        return True, np.zeros(3, dtype=np.float64)
    if n == 1:
        return True, bead_positions[0] + np.array([0, 0, r], dtype=np.float64)

    # Pairwise feasibility: no pair can be farther than 2r apart
    for i in range(n):
        for j in range(i + 1, n):
            delta = bead_positions[j] - bead_positions[i]
            for dim in range(3):
                if pbc[dim] and box_lengths[dim] > 0 and np.isfinite(box_lengths[dim]):
                    delta[dim] -= (
                        np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
                    )
            if np.linalg.norm(delta) > 2.0 * r * (1.0 + tolerance):
                return False, None

    target = _find_equidistant_point(bead_positions, r, pbc, box_lengths)
    distances = np.linalg.norm(target - bead_positions, axis=1)
    max_error = np.max(np.abs(distances - r)) / r
    if max_error <= tolerance:
        return True, target
    return False, None


# =============================================================================
# Overlap Detection (KDTree-accelerated)
# =============================================================================


def _has_overlaps(
    placed_coords, all_coords, excluded_indices, minimum_separation, pbc, box_lengths
):
    """Check if any placed bead overlaps with non-excluded existing beads.

    Parameters
    ----------
    placed_coords : array-like, shape (M, 3)
        Coordinates of beads being placed.
    all_coords : array-like, shape (N, 3)
        All existing bead coordinates.
    excluded_indices : array-like of int or None
        Indices into all_coords to skip (e.g. bonded partners).
    minimum_separation : float
        Overlap threshold distance.
    pbc : array-like of bool, shape (3,)
    box_lengths : array-like, shape (3,)

    Returns
    -------
    bool
    """
    if len(all_coords) == 0 or len(placed_coords) == 0:
        return False

    mask = np.ones(len(all_coords), dtype=bool)
    if excluded_indices is not None:
        for idx in excluded_indices:
            if 0 <= idx < len(all_coords):
                mask[idx] = False

    active_coords = all_coords[mask].astype(np.float64)
    if len(active_coords) == 0:
        return False

    all_pbc = (
        np.all(pbc) and np.all(box_lengths > 0) and np.all(np.isfinite(box_lengths))
    )

    if all_pbc:
        wrapped = active_coords.copy()
        for dim in range(3):
            wrapped[:, dim] %= box_lengths[dim]
        tree = cKDTree(wrapped, boxsize=box_lengths)
        query = np.asarray(placed_coords, dtype=np.float64).copy()
        for dim in range(3):
            query[:, dim] %= box_lengths[dim]
    else:
        tree = cKDTree(active_coords)
        query = np.asarray(placed_coords, dtype=np.float64)

    dists, _ = tree.query(query, k=1)
    return bool(np.any(dists < minimum_separation))


# =============================================================================
# Feasibility Checks
# =============================================================================


def _check_separate_sites_feasibility(
    crosslinker,
    coords,
    group_a,
    group_b,
    crosslink_bond_length,
    tolerance,
    pbc,
    box_lengths,
):
    """Check geometric feasibility for a crosslinker with two distinct physical nodes.

    Verifies that points p0 on sphere(centroid_a, r) and p1 on sphere(centroid_b, r)
    can satisfy |p0 - p1| = internal_distance within tolerance.

    Parameters
    ----------
    crosslinker : CrosslinkerGeometry
    coords : ndarray, shape (N, 3)
        Path coordinates.
    group_a, group_b : list of int
        Backbone bead indices for each physical node.
    crosslink_bond_length : float
    tolerance : float
    pbc : ndarray of bool
    box_lengths : ndarray of float

    Returns
    -------
    bool
    """
    conn_sites = crosslinker.connection_sites
    unique_phys = list(dict.fromkeys(conn_sites))
    if len(unique_phys) < 2:
        return False

    internal_distance = float(
        np.linalg.norm(
            crosslinker.coordinates[unique_phys[0]]
            - crosslinker.coordinates[unique_phys[1]]
        )
    )

    positions_a = np.array([coords[i] for i in group_a], dtype=np.float64)
    positions_b = np.array([coords[i] for i in group_b], dtype=np.float64)

    centroid_a = positions_a.mean(axis=0)
    centroid_b = positions_b.mean(axis=0)

    delta = centroid_b - centroid_a
    for dim in range(3):
        if pbc[dim] and box_lengths[dim] > 0 and np.isfinite(box_lengths[dim]):
            delta[dim] -= np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
    ab_distance = float(np.linalg.norm(delta))

    r = crosslink_bond_length
    min_possible = max(0.0, ab_distance - 2 * r)
    max_possible = ab_distance + 2 * r

    d_min = internal_distance * (1.0 - tolerance)
    d_max = internal_distance * (1.0 + tolerance)

    return min_possible <= d_max and d_min <= max_possible


# =============================================================================
# Candidate Group Finding (KDTree-accelerated)
# =============================================================================


def _find_candidate_groups(
    path,
    crosslinker,
    crosslink_bond_length,
    tolerance,
    excluded_bond_depth,
    max_backbone_degree,
    pbc,
    box_lengths,
    backbone_specs,
):
    """Find groups of backbone nodes that can be connected by the crosslinker.

    Parameters
    ----------
    path : Path
    crosslinker : CrosslinkerGeometry
    crosslink_bond_length : float
    tolerance : float
    excluded_bond_depth : int
    max_backbone_degree : int
    pbc : ndarray of bool
    box_lengths : ndarray of float
    backbone_specs : list or None

    Returns
    -------
    list of list of list of int
        Each candidate is a list (one entry per connection site) of lists
        (bead indices for that site).
    """
    n_conn = len(crosslinker.connection_sites)
    coords = path.coordinates

    bps = _beads_per_site(backbone_specs) if backbone_specs else [1] * n_conn

    # Group connection sites by physical node
    conn_sites = crosslinker.connection_sites
    node_to_conn_indices = {}
    for conn_idx, node in enumerate(conn_sites):
        node_to_conn_indices.setdefault(node, []).append(conn_idx)
    unique_phys_nodes = list(node_to_conn_indices.keys())

    # Gather eligible nodes per connection site
    eligible_per_site = []
    for site_idx in range(n_conn):
        spec = backbone_specs[site_idx] if backbone_specs else None
        eligible = []
        for node in range(len(coords)):
            if not path.bond_graph.has_node(node):
                continue
            if int(path.bond_graph.degree[node]) >= max_backbone_degree:
                continue
            if spec is not None and not _bead_matches_spec(path.beads[node], spec):
                continue
            eligible.append(node)
        eligible_per_site.append(eligible)

    if any(len(e) == 0 for e in eligible_per_site):
        return []

    # Compute search radius
    if len(unique_phys_nodes) == 1:
        search_radius = 2.0 * crosslink_bond_length * (1.0 + tolerance)
    else:
        node_a, node_b = unique_phys_nodes[0], unique_phys_nodes[1]
        internal_dist = float(
            np.linalg.norm(
                crosslinker.coordinates[node_a] - crosslinker.coordinates[node_b]
            )
        )
        search_radius = internal_dist + 2.0 * crosslink_bond_length * (1.0 + tolerance)

    # Dispatch based on topology
    if len(unique_phys_nodes) == 1 and all(b == 1 for b in bps):
        return _pairwise_same_node_search(
            path,
            crosslinker,
            eligible_per_site,
            node_to_conn_indices,
            search_radius,
            crosslink_bond_length,
            tolerance,
            excluded_bond_depth,
            pbc,
            box_lengths,
        )
    elif len(unique_phys_nodes) == 2 and all(b == 1 for b in bps):
        return _pairwise_diff_node_search(
            path,
            crosslinker,
            eligible_per_site,
            node_to_conn_indices,
            unique_phys_nodes,
            search_radius,
            crosslink_bond_length,
            tolerance,
            excluded_bond_depth,
            pbc,
            box_lengths,
        )
    else:
        return _general_candidate_search(
            path,
            crosslinker,
            eligible_per_site,
            node_to_conn_indices,
            unique_phys_nodes,
            bps,
            crosslink_bond_length,
            tolerance,
            excluded_bond_depth,
            pbc,
            box_lengths,
            search_radius,
        )


def _pairwise_same_node_search(
    path,
    crosslinker,
    eligible_per_site,
    node_to_conn_indices,
    search_radius,
    crosslink_bond_length,
    tolerance,
    excluded_bond_depth,
    pbc,
    box_lengths,
):
    """Find candidate pairs when all connections map to the same physical node.

    Each candidate is [[bead_a], [bead_b]] — one bead per connection site.
    """
    coords = path.coordinates
    eligible_set = sorted(set(eligible_per_site[0]) & set(eligible_per_site[1]))
    if len(eligible_set) < 2:
        return []

    eligible_coords = coords[eligible_set].astype(np.float64)
    all_pbc = (
        np.all(pbc) and np.all(box_lengths > 0) and np.all(np.isfinite(box_lengths))
    )

    if all_pbc:
        wrapped = eligible_coords.copy()
        for dim in range(3):
            wrapped[:, dim] %= box_lengths[dim]
        tree = cKDTree(wrapped, boxsize=box_lengths)
    else:
        tree = cKDTree(eligible_coords)

    pairs = tree.query_pairs(search_radius, output_type="ndarray")

    candidates = []
    for p in pairs:
        n0 = eligible_set[p[0]]
        n1 = eligible_set[p[1]]
        if excluded_bond_depth > 0:
            try:
                dist = nx.shortest_path_length(path.bond_graph, n0, n1)
                if dist <= excluded_bond_depth:
                    continue
            except nx.NetworkXNoPath:
                pass
        candidates.append([[n0], [n1]])
        if len(candidates) >= 5000:
            break

    return candidates


def _pairwise_diff_node_search(
    path,
    crosslinker,
    eligible_per_site,
    node_to_conn_indices,
    unique_phys_nodes,
    search_radius,
    crosslink_bond_length,
    tolerance,
    excluded_bond_depth,
    pbc,
    box_lengths,
):
    """Find candidate pairs when connections map to two distinct physical nodes.

    Each candidate is [[bead_a], [bead_b]] — one bead per connection site.
    """
    coords = path.coordinates
    eligible_0 = sorted(set(eligible_per_site[0]))
    eligible_1 = sorted(set(eligible_per_site[1]))

    if not eligible_0 or not eligible_1:
        return []

    coords_0 = coords[eligible_0].astype(np.float64)
    coords_1 = coords[eligible_1].astype(np.float64)

    all_pbc = (
        np.all(pbc) and np.all(box_lengths > 0) and np.all(np.isfinite(box_lengths))
    )

    if all_pbc:
        wrapped_1 = coords_1.copy()
        for dim in range(3):
            wrapped_1[:, dim] %= box_lengths[dim]
        tree_1 = cKDTree(wrapped_1, boxsize=box_lengths)
    else:
        tree_1 = cKDTree(coords_1)

    candidates = []
    for i, n0 in enumerate(eligible_0):
        query_pt = coords_0[i].copy()
        if all_pbc:
            for dim in range(3):
                query_pt[dim] %= box_lengths[dim]

        nearby = tree_1.query_ball_point(query_pt, search_radius)
        for j_local in nearby:
            n1 = eligible_1[j_local]
            if n0 == n1:
                continue
            if excluded_bond_depth > 0:
                try:
                    dist = nx.shortest_path_length(path.bond_graph, n0, n1)
                    if dist <= excluded_bond_depth:
                        continue
                except nx.NetworkXNoPath:
                    pass

            feasible = _check_separate_sites_feasibility(
                crosslinker,
                coords,
                [n0],
                [n1],
                crosslink_bond_length,
                tolerance,
                pbc,
                box_lengths,
            )
            if feasible:
                candidates.append([[n0], [n1]])
                if len(candidates) >= 5000:
                    return candidates

    return candidates


def _general_candidate_search(
    path,
    crosslinker,
    eligible_per_site,
    node_to_conn_indices,
    unique_phys_nodes,
    bps,
    crosslink_bond_length,
    tolerance,
    excluded_bond_depth,
    pbc,
    box_lengths,
    search_radius,
):
    """General candidate search for multi-bead specs or 3+ connection sites.

    Handles cases where connection sites need groups of bonded beads or
    where there are more than two physical connection nodes.

    Parameters
    ----------
    eligible_per_site : list of list of int
        Per connection site, eligible individual bead indices.
    bps : list of int
        Beads-per-site: how many beads each connection site bonds to.
    """
    coords = path.coordinates
    conn_sites = crosslinker.connection_sites
    same_node = len(unique_phys_nodes) == 1

    # Build bonded groups per connection site
    groups_per_site = []
    for site_idx in range(len(conn_sites)):
        n_needed = bps[site_idx]
        eligible = eligible_per_site[site_idx]

        if n_needed == 1:
            groups_per_site.append([[node] for node in eligible])
        else:
            groups = _find_bonded_groups(path, eligible, n_needed)
            groups_per_site.append(groups)

    if any(len(g) == 0 for g in groups_per_site):
        return []

    candidates = []

    if same_node:
        _find_valid_combinations(
            path,
            groups_per_site,
            0,
            [],
            set(),
            crosslink_bond_length,
            tolerance,
            excluded_bond_depth,
            pbc,
            box_lengths,
            coords,
            candidates,
            max_candidates=1000,
        )
    else:
        # Group connection sites by physical node
        phys_node_sites = {}
        for phys_node, conn_indices in node_to_conn_indices.items():
            phys_node_sites[phys_node] = conn_indices

        if len(unique_phys_nodes) == 2:
            node_a, node_b = unique_phys_nodes[0], unique_phys_nodes[1]
            sites_a = phys_node_sites[node_a]
            sites_b = phys_node_sites[node_b]

            combos_a = _get_node_combos(
                path, groups_per_site, sites_a, excluded_bond_depth
            )
            combos_b = _get_node_combos(
                path, groups_per_site, sites_b, excluded_bond_depth
            )

            for combo_a in combos_a:
                beads_a = set()
                for grp in combo_a:
                    beads_a.update(grp)

                for combo_b in combos_b:
                    beads_b = set()
                    for grp in combo_b:
                        beads_b.update(grp)

                    if beads_a & beads_b:
                        continue

                    skip = False
                    if excluded_bond_depth > 0:
                        for ba in beads_a:
                            for bb in beads_b:
                                try:
                                    d = nx.shortest_path_length(path.bond_graph, ba, bb)
                                    if d <= excluded_bond_depth:
                                        skip = True
                                        break
                                except nx.NetworkXNoPath:
                                    pass
                            if skip:
                                break
                    if skip:
                        continue

                    all_beads_a = list(beads_a)
                    all_beads_b = list(beads_b)
                    feasible = _check_separate_sites_feasibility(
                        crosslinker,
                        coords,
                        all_beads_a,
                        all_beads_b,
                        crosslink_bond_length,
                        tolerance,
                        pbc,
                        box_lengths,
                    )
                    if feasible:
                        candidate = [None] * len(conn_sites)
                        for i, si in enumerate(sites_a):
                            candidate[si] = combo_a[i]
                        for i, si in enumerate(sites_b):
                            candidate[si] = combo_b[i]
                        candidates.append(candidate)
                        if len(candidates) >= 1000:
                            return candidates

        elif len(unique_phys_nodes) >= 3:
            node_a, node_b, node_c = (
                unique_phys_nodes[0],
                unique_phys_nodes[1],
                unique_phys_nodes[2],
            )
            sites_a = phys_node_sites[node_a]
            sites_b = phys_node_sites[node_b]
            sites_c = phys_node_sites[node_c]

            combos_a = _get_node_combos(
                path, groups_per_site, sites_a, excluded_bond_depth
            )
            combos_b = _get_node_combos(
                path, groups_per_site, sites_b, excluded_bond_depth
            )
            combos_c = _get_node_combos(
                path, groups_per_site, sites_c, excluded_bond_depth
            )

            for ca in combos_a:
                beads_a = set()
                for grp in ca:
                    beads_a.update(grp)
                for cb in combos_b:
                    beads_b = set()
                    for grp in cb:
                        beads_b.update(grp)
                    if beads_a & beads_b:
                        continue
                    for cc in combos_c:
                        beads_c = set()
                        for grp in cc:
                            beads_c.update(grp)
                        if (beads_a | beads_b) & beads_c:
                            continue

                        skip = False
                        if excluded_bond_depth > 0:
                            all_pairs = [
                                (beads_a, beads_b),
                                (beads_a, beads_c),
                                (beads_b, beads_c),
                            ]
                            for set1, set2 in all_pairs:
                                for b1 in set1:
                                    for b2 in set2:
                                        try:
                                            d = nx.shortest_path_length(
                                                path.bond_graph, b1, b2
                                            )
                                            if d <= excluded_bond_depth:
                                                skip = True
                                                break
                                        except nx.NetworkXNoPath:
                                            pass
                                    if skip:
                                        break
                                if skip:
                                    break
                        if skip:
                            continue

                        candidate = [None] * len(conn_sites)
                        for i, si in enumerate(sites_a):
                            candidate[si] = ca[i]
                        for i, si in enumerate(sites_b):
                            candidate[si] = cb[i]
                        for i, si in enumerate(sites_c):
                            candidate[si] = cc[i]
                        candidates.append(candidate)
                        if len(candidates) >= 1000:
                            return candidates

    return candidates


def _find_bonded_groups(path, eligible, n_needed):
    """Find connected subgraphs of size n_needed among eligible nodes.

    Parameters
    ----------
    path : Path
    eligible : list of int
        Eligible node indices.
    n_needed : int
        Required group size.

    Returns
    -------
    list of list of int
        Each entry is a connected subgraph of the required size.
    """
    eligible_set = set(eligible)
    groups = []
    seen_keys = set()

    if n_needed == 2:
        for node in eligible:
            for neighbor in path.bond_graph.neighbors(node):
                if neighbor in eligible_set and neighbor != node:
                    key = tuple(sorted([node, neighbor]))
                    if key not in seen_keys:
                        seen_keys.add(key)
                        groups.append([node, neighbor])
        return groups

    # General case: DFS from each eligible node
    for start in eligible:
        _dfs_bonded_groups(path, start, eligible_set, n_needed, groups, seen_keys)

    return groups


def _dfs_bonded_groups(path, start, eligible_set, size, groups, seen_keys):
    """DFS to enumerate connected subgraphs of a given size."""
    stack = [([start], {start})]
    while stack:
        current_group, visited = stack.pop()
        if len(current_group) == size:
            key = tuple(sorted(current_group))
            if key not in seen_keys:
                seen_keys.add(key)
                groups.append(list(current_group))
            continue
        for node in current_group:
            for neighbor in path.bond_graph.neighbors(node):
                if neighbor in eligible_set and neighbor not in visited:
                    stack.append((current_group + [neighbor], visited | {neighbor}))


def _find_valid_combinations(
    path,
    groups_per_site,
    site_idx,
    current_selection,
    used_beads,
    crosslink_bond_length,
    tolerance,
    excluded_bond_depth,
    pbc,
    box_lengths,
    coords,
    candidates,
    max_candidates=1000,
):
    """Recursively find valid group combinations across connection sites.

    For same-node crosslinkers: all selected beads must be geometrically
    reachable from a single point at crosslink_bond_length. Groups cannot
    share beads, and excluded_bond_depth is enforced between groups.
    """
    if len(candidates) >= max_candidates:
        return

    if site_idx == len(groups_per_site):
        # All sites assigned — verify geometric feasibility
        all_beads = []
        for grp in current_selection:
            all_beads.extend(grp)

        bead_positions = coords[all_beads].astype(np.float64)

        # Unwrap under PBC relative to first bead
        if len(bead_positions) > 1:
            ref = bead_positions[0].copy()
            for k in range(1, len(bead_positions)):
                delta = bead_positions[k] - ref
                for dim in range(3):
                    if (
                        pbc[dim]
                        and box_lengths[dim] > 0
                        and np.isfinite(box_lengths[dim])
                    ):
                        delta[dim] -= (
                            np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
                        )
                bead_positions[k] = ref + delta

        feasible, _ = _can_reach_all_beads(
            bead_positions, crosslink_bond_length, tolerance, pbc, box_lengths
        )
        if feasible:
            candidates.append([list(grp) for grp in current_selection])
        return

    for group in groups_per_site[site_idx]:
        group_set = set(group)

        if group_set & used_beads:
            continue

        skip = False
        if excluded_bond_depth > 0:
            for bead_new in group:
                for bead_old in used_beads:
                    try:
                        d = nx.shortest_path_length(path.bond_graph, bead_new, bead_old)
                        if d <= excluded_bond_depth:
                            skip = True
                            break
                    except nx.NetworkXNoPath:
                        pass
                if skip:
                    break
        if skip:
            continue

        _find_valid_combinations(
            path,
            groups_per_site,
            site_idx + 1,
            current_selection + [group],
            used_beads | group_set,
            crosslink_bond_length,
            tolerance,
            excluded_bond_depth,
            pbc,
            box_lengths,
            coords,
            candidates,
            max_candidates,
        )


def _get_node_combos(path, groups_per_site, site_indices, excluded_bond_depth):
    """Get valid group combinations for connection sites on one physical node.

    Parameters
    ----------
    path : Path
    groups_per_site : list of list of list of int
    site_indices : list of int
        Which connection site indices belong to this physical node.
    excluded_bond_depth : int

    Returns
    -------
    list of list of list of int
        Each entry is a tuple of groups (one per site_index) with no shared beads.
    """
    if len(site_indices) == 1:
        return [[grp] for grp in groups_per_site[site_indices[0]]]

    combos = []
    _find_node_combos(
        path,
        groups_per_site,
        site_indices,
        0,
        [],
        set(),
        excluded_bond_depth,
        combos,
        max_combos=500,
    )
    return combos


def _find_node_combos(
    path,
    groups_per_site,
    site_indices,
    idx,
    current,
    used_beads,
    excluded_bond_depth,
    combos,
    max_combos,
):
    """Recursively find non-overlapping group combos for sites on one physical node."""
    if len(combos) >= max_combos:
        return

    if idx == len(site_indices):
        combos.append(list(current))
        return

    site = site_indices[idx]
    for group in groups_per_site[site]:
        group_set = set(group)
        if group_set & used_beads:
            continue

        skip = False
        if excluded_bond_depth > 0 and used_beads:
            for bead_new in group:
                for bead_old in used_beads:
                    try:
                        d = nx.shortest_path_length(path.bond_graph, bead_new, bead_old)
                        if d <= excluded_bond_depth:
                            skip = True
                            break
                    except nx.NetworkXNoPath:
                        pass
                if skip:
                    break
        if skip:
            continue

        _find_node_combos(
            path,
            groups_per_site,
            site_indices,
            idx + 1,
            current + [group],
            used_beads | group_set,
            excluded_bond_depth,
            combos,
            max_combos,
        )


# =============================================================================
# Geometric Constraint Solvers
# =============================================================================


def _solve_two_sphere_constraint(A, B, r, d, tolerance):
    """Find points p0 on sphere(A, r) and p1 on sphere(B, r) with |p0 - p1| = d.

    Uses an analytical construction: both points lie in the plane containing
    A and B, offset perpendicular to the AB axis.

    Parameters
    ----------
    A, B : array-like, shape (3,)
        Centers of the two spheres (backbone bead positions).
    r : float
        Sphere radius (crosslink_bond_length).
    d : float
        Required distance between p0 and p1 (internal crosslinker distance).
    tolerance : float
        Fractional tolerance on distances.

    Returns
    -------
    tuple of (p0, p1) or None if infeasible.
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)

    AB = B - A
    ab_dist = np.linalg.norm(AB)

    if ab_dist < 1e-12:
        p0 = A + np.array([0, 0, r])
        p1 = p0 + np.array([d, 0, 0])
        return p0, p1

    # Geometric feasibility
    min_possible = max(0.0, ab_dist - 2 * r)
    max_possible = ab_dist + 2 * r

    if d < min_possible * (1 - tolerance) or d > max_possible * (1 + tolerance):
        return None

    ab_hat = AB / ab_dist

    # Get perpendicular direction
    perp = np.cross(ab_hat, np.array([0, 0, 1], dtype=np.float64))
    if np.linalg.norm(perp) < 1e-10:
        perp = np.cross(ab_hat, np.array([0, 1, 0], dtype=np.float64))
    perp /= np.linalg.norm(perp)

    # Analytical construction:
    # p0 = M - (d/2)*ab_hat + h*perp
    # p1 = M + (d/2)*ab_hat + h*perp
    # where M = (A+B)/2
    # |p0 - A|^2 = (ab_dist/2 - d/2)^2 + h^2 = r^2
    M = (A + B) / 2.0
    half_ab = ab_dist / 2.0
    half_d = d / 2.0

    h_sq = r**2 - (half_ab - half_d) ** 2

    if h_sq < 0:
        return _solve_two_sphere_iterative(A, B, r, d, ab_hat, perp)

    h = np.sqrt(h_sq)

    p0 = M - half_d * ab_hat + h * perp
    p1 = M + half_d * ab_hat + h * perp

    # Verify solution
    d0 = np.linalg.norm(p0 - A)
    d1 = np.linalg.norm(p1 - B)
    dd = np.linalg.norm(p1 - p0)

    if (
        abs(d0 - r) / r <= tolerance
        and abs(d1 - r) / r <= tolerance
        and abs(dd - d) / max(d, 1e-10) <= tolerance
    ):
        return p0, p1

    return _solve_two_sphere_iterative(A, B, r, d, ab_hat, perp)


def _solve_two_sphere_iterative(A, B, r, d, ab_hat, perp):
    """Iterative fallback solver for the two-sphere constraint.

    Searches over a grid of center positions along the AB axis with
    perpendicular offsets, then projects each candidate onto the constraint
    surfaces.

    Parameters
    ----------
    A, B : ndarray, shape (3,)
    r : float
    d : float
    ab_hat : ndarray, shape (3,)
        Unit vector from A to B.
    perp : ndarray, shape (3,)
        Unit vector perpendicular to ab_hat.

    Returns
    -------
    tuple of (p0, p1)
    """
    best_error = np.inf
    best_result = None

    for h_frac in np.linspace(0, 1, 20):
        h = r * h_frac
        for t in np.linspace(0.1, 0.9, 10):
            center = A + t * (B - A) + h * perp
            p0 = center - (d / 2.0) * ab_hat
            p1 = center + (d / 2.0) * ab_hat

            d0 = np.linalg.norm(p0 - A)
            d1 = np.linalg.norm(p1 - B)

            # Project onto spheres
            if d0 > 1e-12:
                p0_adj = A + (p0 - A) * (r / d0)
            else:
                p0_adj = p0
            if d1 > 1e-12:
                p1_adj = B + (p1 - B) * (r / d1)
            else:
                p1_adj = p1

            dd_adj = np.linalg.norm(p1_adj - p0_adj)
            error = abs(dd_adj - d) / max(d, 1e-10)

            if error < best_error:
                best_error = error
                best_result = (p0_adj, p1_adj)

            if error < 0.01:
                return best_result

    if best_result is not None and best_error < 0.5:
        return best_result

    # Last resort: place at bond_length from each, ignore internal distance
    p0 = A + r * perp
    p1 = B + r * perp
    return p0, p1


# =============================================================================
# Rotation Utilities
# =============================================================================


def _rotation_matrix_between_vectors(a, b):
    """Rotation matrix that maps unit vector a to unit vector b.

    Uses Rodrigues' formula via the skew-symmetric cross-product matrix.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    v = np.cross(a, b)
    c = np.dot(a, b)

    if c > 1.0 - 1e-10:
        return np.eye(3, dtype=np.float64)
    if c < -1.0 + 1e-10:
        # 180-degree rotation: find any perpendicular axis
        perp = np.array([1, 0, 0], dtype=np.float64)
        if abs(np.dot(a, perp)) > 0.9:
            perp = np.array([0, 1, 0], dtype=np.float64)
        perp = perp - np.dot(perp, a) * a
        perp /= np.linalg.norm(perp)
        return 2.0 * np.outer(perp, perp) - np.eye(3, dtype=np.float64)

    vx = np.array(
        [
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0],
        ],
        dtype=np.float64,
    )

    return np.eye(3, dtype=np.float64) + vx + vx @ vx / (1.0 + c)


def _rotation_about_axis(axis, angle):
    """Rotation matrix for a given angle about a unit axis (Rodrigues' formula)."""
    axis = np.asarray(axis, dtype=np.float64)
    K = np.array(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ],
        dtype=np.float64,
    )
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def _random_rotation_matrix(rng):
    """Uniformly random rotation matrix via QR decomposition."""
    H = rng.standard_normal((3, 3))
    Q, R = np.linalg.qr(H)
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


# =============================================================================
# Placement with Rotation Sweep
# =============================================================================


def _compute_optimal_placement(
    crosslinker,
    candidate_group,
    path_coords,
    crosslink_bond_length,
    pbc,
    box_lengths,
    tolerance,
):
    """Compute optimal crosslinker coordinates to satisfy bond length constraints.

    Determines target positions for each physical connection node, then
    translates/rotates the crosslinker geometry to align.

    Parameters
    ----------
    crosslinker : CrosslinkerGeometry
    candidate_group : list of list of int
        One list of backbone bead indices per connection site.
    path_coords : ndarray, shape (N, 3)
    crosslink_bond_length : float
    pbc : ndarray of bool
    box_lengths : ndarray of float
    tolerance : float

    Returns
    -------
    placed_coords : ndarray, shape (n_sites, 3) or None if infeasible
    target_node_positions : dict mapping physical node index to position
    """
    conn_sites = crosslinker.connection_sites
    node_to_conn_indices = {}
    for conn_idx, node in enumerate(conn_sites):
        node_to_conn_indices.setdefault(node, []).append(conn_idx)

    unique_phys_nodes = list(node_to_conn_indices.keys())
    target_node_positions = {}

    # Collect all bead indices for PBC unwrapping reference
    all_bead_indices = []
    for group in candidate_group:
        all_bead_indices.extend(group)
    if not all_bead_indices:
        return None, {}

    reference = path_coords[all_bead_indices[0]].astype(np.float64).copy()

    if len(unique_phys_nodes) == 1:
        # Single physical node: all beads must be reachable from one point
        phys_node = unique_phys_nodes[0]

        bead_positions = []
        for group in candidate_group:
            for idx in group:
                pos = path_coords[idx].astype(np.float64)
                delta = pos - reference
                for dim in range(3):
                    if (
                        pbc[dim]
                        and box_lengths[dim] > 0
                        and np.isfinite(box_lengths[dim])
                    ):
                        delta[dim] -= (
                            np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
                        )
                bead_positions.append(reference + delta)
        bead_positions = np.array(bead_positions, dtype=np.float64)

        feasible, target = _can_reach_all_beads(
            bead_positions, crosslink_bond_length, tolerance, pbc, box_lengths
        )
        if not feasible:
            return None, {}
        if target is None:
            target = _find_equidistant_point(
                bead_positions, crosslink_bond_length, pbc, box_lengths
            )

        target_node_positions[phys_node] = target
        node_offset = crosslinker.coordinates[phys_node].astype(np.float64)
        placed_coords = crosslinker.coordinates.astype(np.float64) + (
            target - node_offset
        )

    elif len(unique_phys_nodes) >= 2:
        node_a = unique_phys_nodes[0]
        node_b = unique_phys_nodes[1]
        internal_dist = float(
            np.linalg.norm(
                crosslinker.coordinates[node_a] - crosslinker.coordinates[node_b]
            )
        )

        # Get bead positions per physical node (unwrapped)
        phys_node_beads = {}
        for phys_node, conn_indices in node_to_conn_indices.items():
            positions = []
            for ci in conn_indices:
                if ci < len(candidate_group):
                    for idx in candidate_group[ci]:
                        pos = path_coords[idx].astype(np.float64)
                        delta = pos - reference
                        for dim in range(3):
                            if (
                                pbc[dim]
                                and box_lengths[dim] > 0
                                and np.isfinite(box_lengths[dim])
                            ):
                                delta[dim] -= (
                                    np.round(delta[dim] / box_lengths[dim])
                                    * box_lengths[dim]
                                )
                        positions.append(reference + delta)
            phys_node_beads[phys_node] = (
                np.array(positions, dtype=np.float64) if positions else np.empty((0, 3))
            )

        pos_a = phys_node_beads[node_a]
        pos_b = phys_node_beads[node_b]

        if len(pos_a) == 1 and len(pos_b) == 1:
            result = _solve_two_sphere_constraint(
                pos_a[0], pos_b[0], crosslink_bond_length, internal_dist, tolerance
            )
            if result is None:
                return None, {}
            target_node_positions[node_a] = result[0]
            target_node_positions[node_b] = result[1]
        else:
            # Multi-bead: compute equidistant point per physical node
            for phys_node in unique_phys_nodes:
                bead_pos = phys_node_beads.get(phys_node, np.empty((0, 3)))
                if len(bead_pos) == 0:
                    continue
                feasible, target = _can_reach_all_beads(
                    bead_pos, crosslink_bond_length, tolerance, pbc, box_lengths
                )
                if not feasible:
                    return None, {}
                if target is None:
                    target = _find_equidistant_point(
                        bead_pos, crosslink_bond_length, pbc, box_lengths
                    )
                target_node_positions[phys_node] = target

        # Align crosslinker to target positions
        if len(target_node_positions) == 1:
            node = list(target_node_positions.keys())[0]
            node_offset = crosslinker.coordinates[node].astype(np.float64)
            placed_coords = crosslinker.coordinates.astype(np.float64) + (
                target_node_positions[node] - node_offset
            )
        else:
            src_a = crosslinker.coordinates[node_a].astype(np.float64)
            src_b = crosslinker.coordinates[node_b].astype(np.float64)
            tgt_a = target_node_positions[node_a]
            tgt_b = target_node_positions[node_b]

            src_vec = src_b - src_a
            tgt_vec = tgt_b - tgt_a
            src_len = np.linalg.norm(src_vec)
            tgt_len = np.linalg.norm(tgt_vec)

            if src_len < 1e-12 or tgt_len < 1e-12:
                placed_coords = crosslinker.coordinates.astype(np.float64) + (
                    tgt_a - src_a
                )
            else:
                src_hat = src_vec / src_len
                tgt_hat = tgt_vec / tgt_len
                R = _rotation_matrix_between_vectors(src_hat, tgt_hat)
                centered = crosslinker.coordinates.astype(np.float64) - src_a
                scale = tgt_len / src_len
                placed_coords = (centered * scale) @ R.T + tgt_a

    else:
        return None, {}

    return placed_coords, target_node_positions


def _try_rotated_placements(
    crosslinker,
    candidate_group,
    path_coords,
    crosslink_bond_length,
    pbc,
    box_lengths,
    tolerance,
    minimum_separation,
    excluded_indices,
    n_rotation_samples,
    seed=None,
):
    """Attempt to place crosslinker, trying rotations to avoid overlaps.

    First computes the optimal placement, then if overlaps exist, sweeps
    rotations (axis-aligned for multi-node, random for single-node) to
    find a valid configuration.

    Parameters
    ----------
    crosslinker : CrosslinkerGeometry
    candidate_group : list of list of int
    path_coords : ndarray
    crosslink_bond_length : float
    pbc, box_lengths : ndarray
    tolerance : float
    minimum_separation : float or None
    excluded_indices : list of int
        Bead indices excluded from overlap checks (bonded partners).
    n_rotation_samples : int
    seed : int, optional

    Returns
    -------
    ndarray, shape (n_sites, 3) or None if placement failed.
    """
    base_coords, target_node_positions = _compute_optimal_placement(
        crosslinker,
        candidate_group,
        path_coords,
        crosslink_bond_length,
        pbc,
        box_lengths,
        tolerance,
    )

    if base_coords is None:
        return None

    # Build KDTree for overlap checking
    mask = np.ones(len(path_coords), dtype=bool)
    if excluded_indices is not None:
        for idx in excluded_indices:
            if 0 <= idx < len(path_coords):
                mask[idx] = False
    check_coords = path_coords[mask].astype(np.float64)

    all_pbc = (
        np.all(pbc) and np.all(box_lengths > 0) and np.all(np.isfinite(box_lengths))
    )

    if len(check_coords) > 0:
        if all_pbc:
            wrapped_check = check_coords.copy()
            for dim in range(3):
                wrapped_check[:, dim] %= box_lengths[dim]
            tree = cKDTree(wrapped_check, boxsize=box_lengths)
        else:
            tree = cKDTree(check_coords)
    else:
        return base_coords  # No obstacles

    def _min_clearance(coords):
        query = np.asarray(coords, dtype=np.float64)
        if all_pbc:
            query = query.copy()
            for dim in range(3):
                query[:, dim] %= box_lengths[dim]
        dists, _ = tree.query(query, k=1)
        return float(np.min(dists))

    # Check base placement
    base_clearance = _min_clearance(base_coords)
    if minimum_separation is None or base_clearance >= minimum_separation:
        return base_coords

    # Single-site: vectorized circle sweep
    if crosslinker.n_sites == 1:
        result = _single_site_circle_sweep(
            crosslinker,
            candidate_group,
            path_coords,
            crosslink_bond_length,
            pbc,
            box_lengths,
            n_rotation_samples,
            minimum_separation,
            tree,
            all_pbc,
        )
        if result is not None:
            clearance = _min_clearance(result)
            if clearance >= minimum_separation:
                return result
        return None

    # Multi-site: rotation sweep
    unique_nodes = list(target_node_positions.keys()) if target_node_positions else []

    if len(unique_nodes) >= 2:
        node_a, node_b = unique_nodes[0], unique_nodes[1]
        pos_a = target_node_positions[node_a]
        pos_b = target_node_positions[node_b]
        axis = pos_b - pos_a
        axis_len = np.linalg.norm(axis)
        if axis_len > 1e-10:
            axis_hat = (axis / axis_len).astype(np.float64)
        else:
            axis_hat = np.array([0, 0, 1], dtype=np.float64)
        pivot = ((pos_a + pos_b) / 2.0).astype(np.float64)
        rotation_type = "axis"
    elif len(unique_nodes) == 1:
        pivot = target_node_positions[unique_nodes[0]].astype(np.float64)
        axis_hat = None
        rotation_type = "full"
    else:
        pivot = np.mean(base_coords, axis=0).astype(np.float64)
        axis_hat = None
        rotation_type = "full"

    rng = np.random.default_rng(seed)
    centered = base_coords.astype(np.float64) - pivot

    if rotation_type == "axis":
        angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)
        for angle in angles:
            R = _rotation_about_axis(axis_hat, angle)
            candidate = (centered @ R.T) + pivot
            clearance = _min_clearance(candidate)
            if clearance >= minimum_separation:
                return candidate.astype(np.float32)
    else:
        for _ in range(n_rotation_samples):
            R = _random_rotation_matrix(rng)
            candidate = (centered @ R.T) + pivot
            clearance = _min_clearance(candidate)
            if clearance >= minimum_separation:
                return candidate.astype(np.float32)

    return None


def _single_site_circle_sweep(
    crosslinker,
    candidate_group,
    path_coords,
    crosslink_bond_length,
    pbc,
    box_lengths,
    n_rotation_samples,
    minimum_separation,
    tree,
    all_pbc,
):
    """Vectorized circle sweep for single-site crosslinkers.

    For a single bead equidistant from two backbone beads, the locus of
    valid positions is a circle. This samples that circle and returns the
    point with maximum clearance from existing beads.

    Returns
    -------
    ndarray, shape (1, 3) or None
    """
    bead_positions = []
    for group in candidate_group:
        for idx in group:
            bead_positions.append(path_coords[idx].astype(np.float64))
    bead_positions = np.array(bead_positions, dtype=np.float64)

    if len(bead_positions) < 2:
        target = bead_positions[0] + np.array([0, 0, crosslink_bond_length])
        return target.reshape(1, 3).astype(np.float32)

    # Unwrap under PBC
    ref = bead_positions[0].copy()
    for k in range(1, len(bead_positions)):
        delta = bead_positions[k] - ref
        for dim in range(3):
            if pbc[dim] and box_lengths[dim] > 0 and np.isfinite(box_lengths[dim]):
                delta[dim] -= np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
        bead_positions[k] = ref + delta

    pA, pB = bead_positions[0], bead_positions[1]
    delta = pB - pA
    d = float(np.linalg.norm(delta))
    M = (pA + pB) / 2.0

    if d > 2.0 * crosslink_bond_length:
        return None

    axis = delta / d if d > 1e-10 else np.array([1, 0, 0], dtype=np.float64)
    h = float(np.sqrt(max(0.0, crosslink_bond_length**2 - (d * 0.5) ** 2)))

    if abs(axis[0]) < 0.9:
        perp1 = np.cross(axis, np.array([1, 0, 0], dtype=np.float64))
    else:
        perp1 = np.cross(axis, np.array([0, 1, 0], dtype=np.float64))
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(axis, perp1)
    perp2 /= np.linalg.norm(perp2)

    angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)
    cos_a = np.cos(angles)
    sin_a = np.sin(angles)
    candidates = M[np.newaxis, :] + h * (
        cos_a[:, np.newaxis] * perp1[np.newaxis, :]
        + sin_a[:, np.newaxis] * perp2[np.newaxis, :]
    )

    query_pts = candidates.copy()
    if all_pbc:
        for dim in range(3):
            query_pts[:, dim] %= box_lengths[dim]

    dists, _ = tree.query(query_pts, k=1)

    valid_mask = dists >= minimum_separation
    if np.any(valid_mask):
        valid_indices = np.where(valid_mask)[0]
        best = valid_indices[np.argmax(dists[valid_indices])]
        return candidates[best].reshape(1, 3).astype(np.float32)

    return None


# =============================================================================
# Path Insertion
# =============================================================================


def _insert_crosslinker(path, crosslinker, placed_coords, candidate_group):
    """Insert a placed crosslinker into the path in-place.

    Appends crosslinker beads to path.coordinates and path.beads,
    adds internal crosslinker bonds, and connects each connection site's
    physical node to its assigned backbone beads.

    Parameters
    ----------
    path : Path
    crosslinker : CrosslinkerGeometry
    placed_coords : ndarray, shape (n_sites, 3)
    candidate_group : list of list of int
        One list of backbone bead indices per connection site.
    """
    n_existing = len(path.coordinates)
    n_new = len(placed_coords)

    path.coordinates = np.vstack([path.coordinates, placed_coords.astype(np.float32)])
    path.beads = np.concatenate([path.beads, crosslinker.beads])

    new_node_ids = list(range(n_existing, n_existing + n_new))
    path.bond_graph.add_nodes_from(new_node_ids)

    # Add internal crosslinker bonds
    for u, v in crosslinker.internal_bonds:
        path.bond_graph.add_edge(new_node_ids[u], new_node_ids[v])

    # Add bonds from connection sites to backbone beads
    conn_sites = crosslinker.connection_sites
    for conn_idx, phys_node in enumerate(conn_sites):
        crosslinker_global_id = new_node_ids[phys_node]
        if conn_idx < len(candidate_group):
            for bead_idx in candidate_group[conn_idx]:
                if not path.bond_graph.has_edge(crosslinker_global_id, bead_idx):
                    path.bond_graph.add_edge(crosslinker_global_id, bead_idx)


# =============================================================================
# Main Entry Point
# =============================================================================


def crosslink(
    path,
    crosslinker=None,
    n_crosslinks=1,
    crosslink_bond_length=0.2,
    max_backbone_degree=4,
    tolerance=0.1,
    excluded_bond_depth=2,
    volume_constraint=None,
    minimum_separation=0.1,
    n_rotation_samples=50,
    backbone_name=None,
    seed=None,
    verbose=False,
):
    """Add crosslinks between backbone beads in a path structure.

    Finds eligible backbone bead groups, places crosslinker geometries
    between them at the specified bond length, and modifies the path in-place.

    Parameters
    ----------
    path : Path
        The polymer path to crosslink.
    crosslinker : CrosslinkerGeometry, optional
        The crosslinker geometry to place. Defaults to a single-site
        2-connection crosslinker.
    n_crosslinks : int
        Number of crosslinks to place.
    crosslink_bond_length : float, default 0.2
        Target distance from each crosslinker connection node to its
        backbone bead(s).
    max_backbone_degree : int, default 4
        Maximum graph degree a backbone node may have and remain eligible.
    tolerance : float, default 0.1
        Fractional tolerance on bond lengths (0.1 = ±10%).
    excluded_bond_depth : int, default 2
        Backbone beads within this many bonds of a selected bead are
        excluded from being selected as additional connection points.
    volume_constraint : optional
        Provides PBC information (CuboidConstraint or CylinderConstraint).
    minimum_separation : float, default 0.1
        Minimum distance between placed crosslinker beads and existing beads.
    n_rotation_samples : int, default 50
        Number of rotations to try for overlap avoidance.
    backbone_name : str or list, optional
        Filter eligible backbone beads by name. If a list, must have one
        entry per connection site.
    seed : int, optional
        Random seed for reproducibility.
    verbose : bool, default False
        Print progress information.

    Returns
    -------
    path : Path
        The same path object, modified in-place with crosslinks inserted.

    Raises
    ------
    PathConvergenceError
        If the requested number of crosslinks cannot be placed.
    """
    if crosslinker is None:
        crosslinker = CrosslinkerGeometry.single_site()

    pbc, box_lengths = _get_pbc_info(volume_constraint)
    backbone_specs = (
        _parse_backbone_spec(backbone_name, crosslinker.connection_sites)
        if backbone_name
        else None
    )

    rng = np.random.default_rng(seed)
    n_placed = 0
    max_attempts = n_crosslinks * 100
    attempt = 0

    while n_placed < n_crosslinks and attempt < max_attempts:
        attempt += 1

        candidates = _find_candidate_groups(
            path,
            crosslinker,
            crosslink_bond_length,
            tolerance,
            excluded_bond_depth,
            max_backbone_degree,
            pbc,
            box_lengths,
            backbone_specs,
        )

        if not candidates:
            if verbose:
                print(f"No candidates found after {n_placed} crosslinks.")
            break

        rng.shuffle(candidates)

        placed = False
        for candidate_group in candidates:
            # Beads directly bonded to the crosslinker are excluded from
            # overlap checks (they will be at exactly crosslink_bond_length).
            overlap_excluded = set()
            for group in candidate_group:
                for node in group:
                    overlap_excluded.add(node)
            overlap_excluded = sorted(overlap_excluded)

            placed_coords = _try_rotated_placements(
                crosslinker,
                candidate_group,
                path.coordinates,
                crosslink_bond_length,
                pbc,
                box_lengths,
                tolerance,
                minimum_separation,
                overlap_excluded,
                n_rotation_samples,
                seed=rng.integers(2**31),
            )

            if placed_coords is None:
                continue

            _insert_crosslinker(path, crosslinker, placed_coords, candidate_group)
            n_placed += 1
            placed = True

            if verbose and n_placed % 10 == 0:
                print(f"Placed {n_placed}/{n_crosslinks} crosslinks")
            break

        if not placed:
            break

    if n_placed < n_crosslinks:
        if not candidates:
            raise PathConvergenceError(
                f"Could not find backbone beads matching crosslinker geometry with "
                f"bond_length={crosslink_bond_length} (tolerance=±{tolerance * 100:.0f}%).\n"
                f"max_backbone_degree={max_backbone_degree}, "
                f"excluded_bond_depth={excluded_bond_depth}.\n"
                f"Current crosslinks placed: {n_placed}.\n"
                "Ways to increase crosslinking:\n"
                " - Increase crosslink_bond_length or tolerance\n"
                " - Increase max_backbone_degree\n"
                " - Decrease excluded_bond_depth\n"
                " - Pack at higher density"
            )
        else:
            raise PathConvergenceError(
                f"Could only place {n_placed}/{n_crosslinks} crosslinks after "
                f"{attempt} attempts.\n"
                f"Consider increasing tolerance, n_rotation_samples, or "
                f"decreasing minimum_separation."
            )

    return path
