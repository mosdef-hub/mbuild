"""Crosslinking module for polymer path structures.

Handles placement of crosslinker geometries between backbone beads,
with proper overlap avoidance and bond length enforcement under PBC.
"""

import networkx as nx
import numpy as np

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
        Positions of crosslinker sites (re-centered at origin).
    bond_graph : networkx.Graph, optional
        Internal connectivity.
    bead_name : str or np.ndarray, default "_R"
        Name(s) for each site.
    connection_sites : list of int
        Node indices that form bonds to external beads.
        May contain duplicates (same node bonds to multiple groups).
    """

    def __init__(
        self, coordinates, bond_graph=None, bead_name="_R", connection_sites=None
    ):
        coordinates = np.asarray(coordinates, dtype=np.float32)
        centroid = coordinates.mean(axis=0)
        self.coordinates = coordinates - centroid

        if bond_graph is None:
            bond_graph = nx.Graph()
            for i in range(self.n_sites):
                bond_graph.add_node(i)
        self.bond_graph = bond_graph

        if isinstance(bead_name, str):
            self.beads = np.array([bead_name] * self.n_sites, dtype="U10")
        else:
            self.beads = np.array(bead_name, dtype="U10")

        if connection_sites is None:
            connection_sites = list(range(self.n_sites))
        self.connection_sites = list(connection_sites)

    @property
    def n_sites(self):
        return len(self.coordinates)

    @property
    def n_connections(self):
        return len(self.connection_sites)

    @property
    def unique_connection_sites(self):
        return list(dict.fromkeys(self.connection_sites))

    @property
    def internal_bonds(self):
        return list(self.bond_graph.edges())

    def copy(self):
        """Return a deep copy of this crosslinker geometry."""
        return CrosslinkerGeometry(
            coordinates=self.coordinates.copy(),
            bond_graph=self.bond_graph.copy(),
            bead_name=self.beads.copy(),
            connection_sites=list(self.connection_sites),
        )

    def recenter(self):
        if len(self.coordinates) > 0:
            centroid = np.mean(self.coordinates, axis=0)
            self.coordinates = self.coordinates - centroid

    @classmethod
    def from_path(cls, path, connection_sites):
        """Create from an existing Path.

        Parameters
        ----------
        path : Path
            An existing path defining the structure.
        connection_sites : list of int
            Which nodes bond to external beads.
        """
        return cls(
            coordinates=path.coordinates.copy(),
            bond_graph=path.bond_graph.copy(),
            bead_name=path.beads.copy(),
            connection_sites=connection_sites,
        )

    @classmethod
    def single_site(cls, bead_name="_R", n_connections=2):
        """Single-site (original crosslink behavior)."""
        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        return cls(
            coordinates=coords,
            bead_name=bead_name,
            connection_sites=[0] * n_connections,
        )

    @classmethod
    def equilateral_triangle(
        cls, bond_length=0.27, bead_name="_R", connection_sites=None
    ):
        """Three sites in a flat equilateral triangle with edge = bond_length."""
        R = bond_length / np.sqrt(3)
        angles = [0, 2 * np.pi / 3, 4 * np.pi / 3]
        positions = np.array(
            [[R * np.cos(a), R * np.sin(a), 0.0] for a in angles],
            dtype=np.float32,
        )
        bead_names = (
            np.array([bead_name] * 3, dtype="U10")
            if isinstance(bead_name, str)
            else np.array(bead_name, dtype="U10")
        )
        if connection_sites is None:
            connection_sites = [0, 1, 2]
        G = nx.Graph()
        for i in range(3):
            G.add_node(i)
        G.add_edge(0, 1)
        G.add_edge(1, 2)
        G.add_edge(0, 2)
        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_names,
            connection_sites=connection_sites,
        )

    @classmethod
    def linear(cls, n_sites=2, bond_length=0.27, bead_name="_R", connection_sites=None):
        """Linear chain of sites."""
        positions = np.zeros((n_sites, 3), dtype=np.float32)
        for i in range(n_sites):
            positions[i, 0] = i * bond_length

        bead_names = (
            np.array([bead_name] * n_sites, dtype="U10")
            if isinstance(bead_name, str)
            else np.array(bead_name, dtype="U10")
        )
        if connection_sites is None:
            connection_sites = [0, n_sites - 1]

        G = nx.Graph()
        for i in range(n_sites):
            G.add_node(i)
        for i in range(n_sites - 1):
            G.add_edge(i, i + 1)

        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_names,
            connection_sites=connection_sites,
        )

    @classmethod
    def square(cls, bond_length=0.27, bead_name="_R", connection_sites=None):
        """Four sites in a flat square."""
        half = bond_length / 2.0
        positions = np.array(
            [[half, half, 0], [-half, half, 0], [-half, -half, 0], [half, -half, 0]],
            dtype=np.float32,
        )
        bead_names = (
            np.array([bead_name] * 4, dtype="U10")
            if isinstance(bead_name, str)
            else np.array(bead_name, dtype="U10")
        )
        if connection_sites is None:
            connection_sites = [0, 1, 2, 3]

        G = nx.Graph()
        for i in range(4):
            G.add_node(i)
        for i in range(4):
            G.add_edge(i, (i + 1) % 4)

        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_names,
            connection_sites=connection_sites,
        )

    @classmethod
    def tetrahedral(cls, bond_length=0.27, bead_name="_R", connection_sites=None):
        """Four sites in a regular tetrahedron."""
        positions = np.array(
            [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=np.float32
        )
        current_dist = np.linalg.norm(positions[0] - positions[1])
        positions *= bond_length / current_dist

        bead_names = (
            np.array([bead_name] * 4, dtype="U10")
            if isinstance(bead_name, str)
            else np.array(bead_name, dtype="U10")
        )
        if connection_sites is None:
            connection_sites = [0, 1, 2, 3]

        G = nx.Graph()
        for i in range(4):
            G.add_node(i)
        for i in range(4):
            for j in range(i + 1, 4):
                G.add_edge(i, j)

        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_names,
            connection_sites=connection_sites,
        )

    @classmethod
    def trigonal_bipyramidal(
        cls, bond_length=0.27, bead_name="_R", connection_sites=None
    ):
        """Five sites in trigonal bipyramidal arrangement (0-2 equatorial, 3-4 axial)."""
        R_eq = bond_length / np.sqrt(3)
        eq_angles = [0, 2 * np.pi / 3, 4 * np.pi / 3]
        equatorial = [[R_eq * np.cos(a), R_eq * np.sin(a), 0.0] for a in eq_angles]
        axial_dist = bond_length * np.sqrt(2.0 / 3.0)
        axial = [[0.0, 0.0, axial_dist], [0.0, 0.0, -axial_dist]]
        positions = np.array(equatorial + axial, dtype=np.float32)

        bead_names = (
            np.array([bead_name] * 5, dtype="U10")
            if isinstance(bead_name, str)
            else np.array(bead_name, dtype="U10")
        )
        if connection_sites is None:
            connection_sites = [0, 1, 2, 3, 4]

        G = nx.Graph()
        for i in range(5):
            G.add_node(i)
        G.add_edge(0, 1)
        G.add_edge(1, 2)
        G.add_edge(0, 2)
        for eq in range(3):
            G.add_edge(eq, 3)
            G.add_edge(eq, 4)

        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_names,
            connection_sites=connection_sites,
        )

    @classmethod
    def pentagon(cls, bond_length=0.27, bead_name="_R", connection_sites=None):
        """Five sites in a flat regular pentagon."""
        R = bond_length / (2 * np.sin(np.pi / 5))
        angles = [2 * np.pi * i / 5 for i in range(5)]
        positions = np.array(
            [[R * np.cos(a), R * np.sin(a), 0.0] for a in angles], dtype=np.float32
        )

        bead_names = (
            np.array([bead_name] * 5, dtype="U10")
            if isinstance(bead_name, str)
            else np.array(bead_name, dtype="U10")
        )
        if connection_sites is None:
            connection_sites = [0, 1, 2, 3, 4]

        G = nx.Graph()
        for i in range(5):
            G.add_node(i)
        for i in range(5):
            G.add_edge(i, (i + 1) % 5)

        return cls(
            coordinates=positions,
            bond_graph=G,
            bead_name=bead_names,
            connection_sites=connection_sites,
        )

    @classmethod
    def from_edges(cls, coordinates, edges, bead_name="_R", connection_sites=None):
        """Create from explicit positions and edge list.

        Parameters
        ----------
        coordinates : array-like, shape (N, 3)
        edges : list of tuple(int, int)
        bead_name : str or list of str
        connection_sites : list of int, optional
        """
        coordinates = np.asarray(coordinates, dtype=np.float32)
        n = len(coordinates)
        bead_names = (
            np.array([bead_name] * n, dtype="U10")
            if isinstance(bead_name, str)
            else np.array(bead_name, dtype="U10")
        )
        if connection_sites is None:
            connection_sites = list(range(n))

        G = nx.Graph()
        for i in range(n):
            G.add_node(i)
        for i, j in edges:
            G.add_edge(i, j)

        return cls(
            coordinates=coordinates,
            bond_graph=G,
            bead_name=bead_names,
            connection_sites=connection_sites,
        )


# =============================================================================
# PBC Utilities
# =============================================================================


def _pbc_delta(a, b, box_lengths, pbc):
    """Minimum-image displacement: a - b."""
    delta = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    for dim in range(3):
        if pbc[dim] and box_lengths[dim] > 0:
            delta[..., dim] -= (
                np.round(delta[..., dim] / box_lengths[dim]) * box_lengths[dim]
            )
    return delta.astype(np.float32)


def _pbc_distance(a, b, box_lengths, pbc):
    """Minimum-image scalar distance."""
    delta = _pbc_delta(a, b, box_lengths, pbc)
    return np.linalg.norm(delta, axis=-1)


def _wrap_positions(positions, box_lengths, pbc):
    """Wrap positions into periodic box centered at origin."""
    positions = positions.copy()
    for dim in range(3):
        if pbc[dim] and box_lengths[dim] > 0:
            half = box_lengths[dim] / 2
            positions[:, dim] = ((positions[:, dim] + half) % box_lengths[dim]) - half
    return positions


# =============================================================================
# Backbone Specification Parsing
# =============================================================================


def _parse_backbone_spec(backbone_name, connection_sites):
    """Parse backbone_name into per-connection-site specifications.

    Parameters
    ----------
    backbone_name : str or tuple
        - String: every connection site bonds to one bead of that type.
        - Tuple of length == len(connection_sites): per-site spec.
            - String entry: one backbone bead.
            - Tuple entry ("_B", "_B"): multiple neighboring backbone beads.

    Returns
    -------
    specs : list of list[str]
        specs[i] is a list of bead names that connection site i bonds to.
    """
    n_conn = len(connection_sites)

    if isinstance(backbone_name, str):
        return [[backbone_name]] * n_conn

    if not isinstance(backbone_name, (tuple, list)):
        raise ValueError(
            f"backbone_name must be str or tuple, got {type(backbone_name)}"
        )

    if len(backbone_name) != n_conn:
        raise ValueError(
            f"backbone_name length ({len(backbone_name)}) must match "
            f"number of connection sites ({n_conn})."
        )

    specs = []
    for entry in backbone_name:
        if isinstance(entry, str):
            specs.append([entry])
        elif isinstance(entry, (tuple, list)):
            specs.append(list(entry))
        else:
            raise ValueError(f"Entry must be str or tuple, got {type(entry)}")

    return specs


# =============================================================================
# Geometry & Graph Helpers
# =============================================================================


def _get_excluded_indices(path, node_indices, excluded_bond_depth):
    """Get indices within excluded_bond_depth bonds of any node in node_indices."""
    excluded = set(node_indices)
    frontier = set(node_indices)
    for _ in range(excluded_bond_depth):
        next_frontier = set()
        for node in frontier:
            if node in path.bond_graph:
                for neighbor in path.bond_graph.neighbors(node):
                    if neighbor not in excluded:
                        excluded.add(neighbor)
                        next_frontier.add(neighbor)
        frontier = next_frontier
    return excluded


def _get_perpendicular(v):
    """Return a unit vector perpendicular to v."""
    v = np.asarray(v, dtype=np.float64)
    if abs(v[0]) < 0.9:
        perp = np.cross(v, [1, 0, 0])
    else:
        perp = np.cross(v, [0, 1, 0])
    norm = np.linalg.norm(perp)
    return (perp / max(norm, 1e-10)).astype(np.float32)


# =============================================================================
# Candidate Search
# =============================================================================


def _find_neighbor_groups(path, target_names, max_backbone_degree):
    """Find bonded groups of backbone beads matching target_names spec.

    For ("_B", "_B"), finds pairs of directly bonded beads where both
    have degree < max_backbone_degree and names match.

    Returns
    -------
    groups : list of list[int]
    """
    n_needed = len(target_names)

    all_eligible = set()
    for i in range(len(path.beads)):
        if path.bond_graph.degree(i) < max_backbone_degree:
            if path.beads[i] in target_names:
                all_eligible.add(i)

    if n_needed == 1:
        return [[i] for i in all_eligible if path.beads[i] == target_names[0]]

    if n_needed == 2:
        groups = []
        seen = set()
        for i in all_eligible:
            for j in path.bond_graph.neighbors(i):
                if j not in all_eligible:
                    continue
                key = (min(i, j), max(i, j))
                if key in seen:
                    continue
                seen.add(key)
                if (
                    path.beads[i] == target_names[0]
                    and path.beads[j] == target_names[1]
                ):
                    groups.append([i, j])
                elif (
                    path.beads[i] == target_names[1]
                    and path.beads[j] == target_names[0]
                ):
                    groups.append([j, i])
        return groups

    # General case for n_needed > 2
    groups = []
    found_keys = set()
    for start in all_eligible:
        _dfs_groups(
            path, start, all_eligible, target_names, n_needed, groups, found_keys
        )
    return groups


def _dfs_groups(path, start, eligible, target_names, size, groups, found_keys):
    """DFS to find connected groups of given size from start node."""
    stack = [([start], {start})]
    while stack:
        current_group, visited = stack.pop()
        if len(current_group) == size:
            key = tuple(sorted(current_group))
            if key not in found_keys:
                group_names = sorted(path.beads[n] for n in current_group)
                if group_names == sorted(target_names):
                    found_keys.add(key)
                    groups.append(list(current_group))
            continue
        last = current_group[-1]
        for neighbor in path.bond_graph.neighbors(last):
            if neighbor in eligible and neighbor not in visited:
                stack.append((current_group + [neighbor], visited | {neighbor}))


def _all_beads_from_groups(candidate_group):
    """Flatten a candidate group into a single list of all bead indices."""
    all_beads = []
    for group in candidate_group:
        all_beads.extend(group)
    return all_beads


def _can_reach_all_beads(
    bead_positions, crosslink_bond_length, tolerance, pbc, box_lengths
):
    """Check if a point exists within tolerance of crosslink_bond_length from all beads."""
    n = len(bead_positions)
    bead_positions = np.asarray(bead_positions, dtype=np.float64)

    if n == 0:
        return True, np.zeros(3, dtype=np.float32)

    if n == 1:
        offset = np.array([0.0, 0.0, crosslink_bond_length], dtype=np.float64)
        return True, (bead_positions[0] + offset).astype(np.float32)

    centroid = bead_positions.mean(axis=0)
    spread = bead_positions - centroid
    principal = (
        spread[0] if np.linalg.norm(spread[0]) > 1e-10 else np.array([1.0, 0.0, 0.0])
    )
    perp = _get_perpendicular(principal).astype(np.float64)

    # Multiple starting points — centroid FIRST (finds minimum-radius in-plane solution)
    starts = [
        centroid.copy(),
        centroid + perp * crosslink_bond_length,
        centroid - perp * crosslink_bond_length,
        centroid + perp * crosslink_bond_length * 0.1,
    ]

    best_point = None
    best_radius_error = float("inf")

    for start in starts:
        point = start.copy()

        # Equidistant solver: converge to point with equal distances to all beads
        for iteration in range(300):
            distances = np.array(
                [np.linalg.norm(point - bead_positions[i]) for i in range(n)]
            )
            mean_dist = distances.mean()
            if mean_dist < 1e-12:
                point += perp * 0.01
                continue
            shift = np.zeros(3, dtype=np.float64)
            for i in range(n):
                delta = point - bead_positions[i]
                dist = distances[i]
                if dist < 1e-12:
                    delta = perp * 0.001
                    dist = np.linalg.norm(delta)
                desired = bead_positions[i] + delta * (mean_dist / dist)
                shift += desired - point
            shift /= n
            point += shift
            if np.linalg.norm(shift) < 1e-10:
                break

        # Check quality
        distances = np.array(
            [np.linalg.norm(point - bead_positions[i]) for i in range(n)]
        )
        natural_radius = distances.mean()
        max_spread = np.max(np.abs(distances - natural_radius))

        if max_spread / max(natural_radius, 1e-10) > 0.05:
            continue  # didn't converge to equidistant point

        radius_error = (
            abs(natural_radius - crosslink_bond_length) / crosslink_bond_length
        )
        if radius_error < best_radius_error:
            best_radius_error = radius_error
            best_point = point.copy()

    # Also try fixed-radius solver from centroid
    for start in starts[:2]:
        point = start.copy()
        for iteration in range(300):
            shift = np.zeros(3, dtype=np.float64)
            for i in range(n):
                delta = point - bead_positions[i]
                dist = np.linalg.norm(delta)
                if dist < 1e-12:
                    delta = perp * 0.001
                    dist = np.linalg.norm(delta)
                desired = bead_positions[i] + delta * (crosslink_bond_length / dist)
                shift += desired - point
            shift /= n
            point += shift
            if np.linalg.norm(shift) < crosslink_bond_length * 1e-8:
                break

        distances = np.array(
            [np.linalg.norm(point - bead_positions[i]) for i in range(n)]
        )
        max_err = (
            np.max(np.abs(distances - crosslink_bond_length)) / crosslink_bond_length
        )
        if max_err < best_radius_error:
            best_radius_error = max_err
            best_point = point.copy()

    if best_point is not None and best_radius_error <= tolerance:
        return True, best_point.astype(np.float32)

    return False, None


def _find_candidate_groups(
    path,
    crosslinker,
    backbone_specs,
    crosslink_bond_length,
    tolerance,
    excluded_bond_depth,
    max_backbone_degree,
    pbc,
    box_lengths,
    rng,
):
    """Find groups of backbone nodes that can be connected by the crosslinker.

    Each candidate is a list of groups (one per connection site), where each
    group is a list of backbone node indices.

    The key constraint: for each connection site, the crosslinker node at that
    site must be placeable at crosslink_bond_length from EACH bead in its group.
    """
    n_conn = len(crosslinker.connection_sites)
    coords = path.coordinates

    # --- Gather eligible groups per connection site ---
    eligible_per_site = []
    for site_idx in range(n_conn):
        spec = backbone_specs[site_idx]
        groups = _find_neighbor_groups(path, spec, max_backbone_degree)
        eligible_per_site.append(groups)

    if any(len(e) == 0 for e in eligible_per_site):
        return []

    for e in eligible_per_site:
        rng.shuffle(e)

    # --- Determine which connection sites share the same physical node ---
    # Groups that share a connection site node have coupled constraints
    conn_site_nodes = crosslinker.connection_sites  # e.g., [0, 0] or [0, 1]

    # Find unique physical nodes and which connection indices map to them
    node_to_conn_indices = {}
    for conn_idx, node in enumerate(conn_site_nodes):
        node_to_conn_indices.setdefault(node, []).append(conn_idx)

    # --- Search for valid combinations ---
    if n_conn == 2:
        return _pairwise_candidate_search(
            path,
            crosslinker,
            eligible_per_site,
            node_to_conn_indices,
            crosslink_bond_length,
            tolerance,
            excluded_bond_depth,
            pbc,
            box_lengths,
            rng,
        )

    return _general_candidate_search(
        path,
        crosslinker,
        eligible_per_site,
        node_to_conn_indices,
        crosslink_bond_length,
        tolerance,
        excluded_bond_depth,
        pbc,
        box_lengths,
        rng,
    )


def _pairwise_candidate_search(
    path,
    crosslinker,
    eligible_per_site,
    node_to_conn_indices,
    crosslink_bond_length,
    tolerance,
    excluded_bond_depth,
    pbc,
    box_lengths,
    rng,
):
    """Search for valid 2-connection-site candidates.

    Handles two cases:
    1. Both connection sites are DIFFERENT physical nodes → check that each
       group is reachable from its own connection site, plus geometric
       compatibility between sites.
    2. Both connection sites are the SAME physical node → all beads from
       both groups must be reachable from a single point.
    """
    coords = path.coordinates
    candidates = []

    # Are both connection sites the same physical node?
    conn_sites = crosslinker.connection_sites
    same_node = conn_sites[0] == conn_sites[1]

    for group_a in eligible_per_site[0]:
        excluded_a = _get_excluded_indices(path, group_a, excluded_bond_depth)

        for group_b in eligible_per_site[1]:
            # No overlap between groups
            if set(group_b) & set(group_a):
                continue
            # Exclusion zone
            if set(group_b) & excluded_a:
                continue
            if set(group_a) & _get_excluded_indices(path, group_b, excluded_bond_depth):
                continue

            # --- Geometric feasibility check ---
            if same_node:
                # Both groups bond to the SAME crosslinker bead.
                # All beads from both groups must be reachable from one point.
                all_bead_indices = list(group_a) + list(group_b)
                all_bead_positions = coords[all_bead_indices].copy()
                feasible, _ = _can_reach_all_beads(
                    all_bead_positions,
                    crosslink_bond_length,
                    tolerance,
                    pbc,
                    box_lengths,
                )
            else:
                # Different physical nodes. Each group independently checked,
                # plus pairwise distance between connection sites must match geometry.
                feasible = _check_separate_sites_feasibility(
                    crosslinker,
                    coords,
                    group_a,
                    group_b,
                    crosslink_bond_length,
                    tolerance,
                    pbc,
                    box_lengths,
                )

            if feasible:
                candidates.append([group_a, group_b])

            if len(candidates) >= 200:
                return candidates

    return candidates


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
    """Check feasibility when connection sites are different physical nodes.

    Verifies:
    1. Each group's beads can be reached from a point at crosslink_bond_length.
    2. The distance between the two reach-points is compatible with the
       rigid crosslinker geometry (distance between connection site nodes).
    """
    conn_sites = crosslinker.connection_sites
    conn_positions = crosslinker.coordinates[conn_sites]
    internal_distance = float(np.linalg.norm(conn_positions[0] - conn_positions[1]))

    # Check group_a reachability
    positions_a = coords[group_a].copy()
    ok_a, point_a = _can_reach_all_beads(
        positions_a, crosslink_bond_length, tolerance, pbc, box_lengths
    )
    if not ok_a:
        return False

    # Check group_b reachability
    positions_b = coords[group_b].copy()
    ok_b, point_b = _can_reach_all_beads(
        positions_b, crosslink_bond_length, tolerance, pbc, box_lengths
    )
    if not ok_b:
        return False

    # Check geometric feasibility via triangle inequality.
    # A rigid rod of length d (internal distance) must have one end on
    # sphere(centroid_a, r) and the other on sphere(centroid_b, r).
    # This is feasible when: (D-d)/2 <= r, i.e., h² = r² - ((D-d)/2)² >= 0
    centroid_a = np.mean(coords[group_a], axis=0).astype(np.float64)
    centroid_b = np.mean(coords[group_b], axis=0).astype(np.float64)
    D = float(_pbc_distance(centroid_a, centroid_b, box_lengths, pbc))
    r = crosslink_bond_length * (1 + tolerance)

    half_diff = abs(D - internal_distance) / 2.0
    if half_diff > r:
        return False

    return True


def _general_candidate_search(
    path,
    crosslinker,
    eligible_per_site,
    node_to_conn_indices,
    crosslink_bond_length,
    tolerance,
    excluded_bond_depth,
    pbc,
    box_lengths,
    rng,
):
    """General search for n-connection-site crosslinkers."""
    coords = path.coordinates
    candidates = []
    max_tries = 1000
    n_conn = len(crosslinker.connection_sites)

    for _ in range(max_tries):
        idx_0 = rng.integers(len(eligible_per_site[0]))
        group_0 = eligible_per_site[0][idx_0]
        excluded = _get_excluded_indices(path, group_0, excluded_bond_depth)
        selected = [group_0]

        # Track accumulated beads per physical node
        node_to_beads = {}  # phys_node_id -> list of bead indices
        phys_0 = crosslinker.connection_sites[0]
        node_to_beads[phys_0] = list(group_0)

        all_nodes = set(group_0)
        valid = True

        for site_idx in range(1, n_conn):
            phys_node = crosslinker.connection_sites[site_idx]
            found = False

            iterGroups = list(rng.permutation(eligible_per_site[site_idx]))
            for group in iterGroups:
                if set(group) & (all_nodes | excluded):
                    continue

                # --- Accumulate ALL beads for this physical node ---
                accumulated = node_to_beads.get(phys_node, []) + list(group)
                accumulated_positions = coords[accumulated].copy()

                # --- Check if ALL accumulated beads are reachable from ONE point ---
                feasible, _ = _can_reach_all_beads(
                    accumulated_positions,
                    crosslink_bond_length,
                    tolerance,
                    pbc,
                    box_lengths,
                )
                if not feasible:
                    continue

                # --- Passed: update state ---
                selected.append(group)
                all_nodes.update(group)
                node_to_beads.setdefault(phys_node, []).extend(group)
                excluded |= _get_excluded_indices(path, group, excluded_bond_depth)
                found = True
                break

            if not found:
                valid = False
                # print(f"Bad Candidate: {idx_0=}\t{group_0=}\t{selected=}")
                break

        if valid:
            candidates.append(selected)
            # print(f"Good Candidate: {idx_0=}\t{group_0=}\t{selected=}")
            if len(candidates) >= 50:
                break

    # --- Deduplicate candidates ---
    # Group indices by physical node, sort within each group, then use as hashable key
    unique_candidates = []
    seen = set()
    for candidate in candidates:
        # Build normalized key: for each physical node, collect and sort all bead indices
        node_to_sorted_beads = {}
        for site_idx, group in enumerate(candidate):
            phys_node = crosslinker.connection_sites[site_idx]
            node_to_sorted_beads.setdefault(phys_node, []).extend(int(x) for x in group)
        # Sort beads within each physical node, then sort by physical node
        key = tuple(
            tuple(sorted(beads)) for _, beads in sorted(node_to_sorted_beads.items())
        )
        if key not in seen:
            seen.add(key)
            unique_candidates.append(candidate)

    return unique_candidates


# =============================================================================
# Crosslinker Placement
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
    """Compute optimal crosslinker placement satisfying bond length constraints.

    For each PHYSICAL node in the crosslinker, collects all backbone beads
    that bond to it (across possibly multiple connection sites) and finds
    the position that minimizes bond length error.

    Returns
    -------
    placed_coords : np.ndarray, shape (n_sites, 3)
        Crosslinker coordinates after optimal placement.
    target_node_positions : dict of {int: np.ndarray}
        Computed position for each physical crosslinker node that has
        external bonds. Keys are node indices in the crosslinker.
    """
    conn_sites = crosslinker.connection_sites

    # --- Group backbone beads by physical crosslinker node ---
    # node_to_beads[phys_node] = list of all backbone bead indices bonding to it
    node_to_beads = {}
    for conn_idx, phys_node in enumerate(conn_sites):
        node_to_beads.setdefault(phys_node, [])
        node_to_beads[phys_node].extend(candidate_group[conn_idx])

    # --- For each physical node, find optimal position ---
    # Use first group's first bead as PBC reference
    reference = path_coords[candidate_group[0][0]].copy()

    target_node_positions = {}
    for phys_node, bead_indices in node_to_beads.items():
        # Unwrap bead positions relative to reference
        bead_positions = []
        for idx in bead_indices:
            delta = _pbc_delta(path_coords[idx], reference, box_lengths, pbc)
            bead_positions.append(reference + delta)
        bead_positions = np.array(bead_positions, dtype=np.float64)

        # Find point equidistant (at crosslink_bond_length) from all beads
        # If impossible, find natural equidistant point within tolerance
        _, target = _can_reach_all_beads(
            bead_positions, crosslink_bond_length, tolerance, pbc, box_lengths
        )
        if target is None:
            target = _find_equidistant_point(
                bead_positions, crosslink_bond_length, pbc, box_lengths
            )
        target_node_positions[phys_node] = target

    # --- Place the crosslinker based on computed node positions ---
    unique_nodes = list(target_node_positions.keys())

    if len(unique_nodes) == 1:
        # Single physical node: just translate the crosslinker
        node = unique_nodes[0]
        node_offset = crosslinker.coordinates[node]
        centroid_shift = target_node_positions[node] - node_offset
        placed_coords = crosslinker.coordinates.astype(np.float64) + centroid_shift
    elif len(unique_nodes) >= 2:
        # Get unwrapped bead centroids for each physical node
        node_bead_centroids = {}
        for phys_node, bead_indices in node_to_beads.items():
            bead_positions = []
            for idx in bead_indices:
                delta = _pbc_delta(path_coords[idx], reference, box_lengths, pbc)
                bead_positions.append(reference + delta)
            node_bead_centroids[phys_node] = np.mean(bead_positions, axis=0).astype(
                np.float64
            )

        if len(unique_nodes) == 2:
            n0, n1 = unique_nodes[0], unique_nodes[1]
            bead_a = node_bead_centroids[n0]
            bead_b = node_bead_centroids[n1]

            # Internal distance between connection sites (rigid)
            conn_0_local = crosslinker.coordinates[n0].astype(np.float64)
            conn_1_local = crosslinker.coordinates[n1].astype(np.float64)
            d = float(np.linalg.norm(conn_1_local - conn_0_local))

            # Backbone distance
            D = float(np.linalg.norm(bead_b - bead_a))
            r = float(crosslink_bond_length)

            # Geometric solution: midpoint of targets is at midpoint of beads,
            # offset perpendicular by h = sqrt(r² - ((D-d)/2)²)
            half_diff = (D - d) / 2.0
            h_sq = r * r - half_diff * half_diff
            h = np.sqrt(max(h_sq, 0.0))

            # Coordinate frame
            midpoint_bb = (bead_a + bead_b) / 2.0
            bb_axis = (bead_b - bead_a) / max(D, 1e-10)
            perp = _get_perpendicular(bb_axis).astype(np.float64)

            # Target positions for the two connection sites
            target_0 = midpoint_bb - (d / 2.0) * bb_axis + h * perp
            target_1 = midpoint_bb + (d / 2.0) * bb_axis + h * perp

            target_node_positions[n0] = target_0.astype(np.float32)
            target_node_positions[n1] = target_1.astype(np.float32)

        else:
            # 3+ nodes: use independently computed targets (best effort)
            pass  # target_node_positions already computed above

        # Kabsch alignment: map crosslinker node positions → target positions
        source_points = np.array(
            [crosslinker.coordinates[n].astype(np.float64) for n in unique_nodes]
        )
        target_points = np.array(
            [target_node_positions[n].astype(np.float64) for n in unique_nodes]
        )

        src_centroid = source_points.mean(axis=0)
        tgt_centroid = target_points.mean(axis=0)

        src_centered = source_points - src_centroid
        tgt_centered = target_points - tgt_centroid

        H = src_centered.T @ tgt_centered
        U, S, Vt = np.linalg.svd(H)
        det = np.linalg.det(Vt.T @ U.T)
        sign_matrix = np.diag([1.0, 1.0, np.sign(det) if abs(det) > 1e-10 else 1.0])
        R = Vt.T @ sign_matrix @ U.T

        all_centered = crosslinker.coordinates.astype(np.float64) - src_centroid
        placed_coords = (all_centered @ R.T) + tgt_centroid

        # Update target_node_positions to actual placed positions
        for n in unique_nodes:
            target_node_positions[n] = placed_coords[n].astype(np.float32)
    else:
        # No connections (shouldn't happen, but be safe)
        placed_coords = crosslinker.coordinates.astype(np.float64)

    return placed_coords.astype(np.float32), target_node_positions


def _find_equidistant_point(bead_positions, target_distance, pbc, box_lengths):
    """Find a point at target_distance from all given bead positions.

    Uses iterative projection onto the intersection of spheres of radius
    target_distance centered at each bead.

    Parameters
    ----------
    bead_positions : np.ndarray, shape (N, 3)
        Positions of backbone beads (already PBC-unwrapped).
    target_distance : float
        Desired distance from the result to each bead.

    Returns
    -------
    point : np.ndarray, shape (3,)
    """
    n = len(bead_positions)
    bead_positions = np.asarray(bead_positions, dtype=np.float64)

    if n == 1:
        offset = np.array([0.0, 0.0, target_distance], dtype=np.float64)
        return (bead_positions[0] + offset).astype(np.float32)

    centroid = bead_positions.mean(axis=0)
    spread = bead_positions - centroid
    principal = (
        spread[0] if np.linalg.norm(spread[0]) > 1e-10 else np.array([1.0, 0.0, 0.0])
    )
    perp = _get_perpendicular(principal).astype(np.float64)
    point = centroid + perp * target_distance

    for iteration in range(100):
        # Compute gradient: move toward the surface of each sphere
        shift = np.zeros(3, dtype=np.float64)
        for i in range(n):
            delta = point - bead_positions[i]
            dist = np.linalg.norm(delta)
            if dist < 1e-12:
                # Point is at bead position; nudge it
                delta = np.array([1e-6, 1e-6, 1e-6], dtype=np.float64)
                dist = np.linalg.norm(delta)
            # Project point onto sphere i: move to radius target_distance
            desired = bead_positions[i] + delta * (target_distance / dist)
            shift += desired - point

        shift /= n
        point += shift

        if np.linalg.norm(shift) < target_distance * 1e-8:
            break

    return point.astype(np.float32)


# =============================================================================
# Rotation Matrix Utilities
# =============================================================================


def _kabsch_rotation(P, Q):
    """Optimal rotation aligning source P to target Q (both centered).

    Returns R such that R @ P ≈ Q in least-squares sense.
    """
    P = np.asarray(P, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, np.sign(d) if abs(d) > 1e-10 else 1.0])
    R = Vt.T @ sign_matrix @ U.T
    return R


def _rotation_between_vectors(src, tgt):
    """Rotation matrix taking unit vector src to unit vector tgt."""
    src = np.asarray(src, dtype=np.float64)
    tgt = np.asarray(tgt, dtype=np.float64)
    v = np.cross(src, tgt)
    c = float(np.dot(src, tgt))
    if c < -0.9999:
        perp = _get_perpendicular(src)
        return (2 * np.outer(perp, perp) - np.eye(3)).astype(np.float32)
    s = np.linalg.norm(v)
    if s < 1e-10:
        return np.eye(3, dtype=np.float32)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R = np.eye(3) + vx + (vx @ vx) * (1.0 / (1.0 + c))
    return R.astype(np.float32)


def _rotate_around_axis(points, axis, angle):
    """Rodrigues rotation of points around axis through the origin."""
    axis = np.asarray(axis, dtype=np.float64)
    norm = np.linalg.norm(axis)
    if norm < 1e-10:
        return points.copy()
    axis = axis / norm
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    K = np.array(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ]
    )
    R = np.eye(3) * cos_a + sin_a * K + (1 - cos_a) * np.outer(axis, axis)
    return (R @ np.asarray(points, dtype=np.float64).T).T.astype(np.float32)


# =============================================================================
# Overlap Checking
# =============================================================================


def _has_overlaps(
    placed_coords, all_coords, excluded_indices, minimum_separation, pbc, box_lengths
):
    """Check if any placed bead overlaps with non-excluded existing beads."""
    if len(all_coords) == 0:
        return False

    mask = np.ones(len(all_coords), dtype=bool)
    for idx in excluded_indices:
        if 0 <= idx < len(all_coords):
            mask[idx] = False
    active_coords = all_coords[mask]

    if len(active_coords) == 0:
        return False

    for point in placed_coords:
        deltas = _pbc_delta(point[np.newaxis, :], active_coords, box_lengths, pbc)
        distances_sq = np.sum(deltas**2, axis=1)
        if np.any(distances_sq < minimum_separation**2):
            return True

    return False


# =============================================================================
# Bond Length Validation
# =============================================================================


def _compute_max_bond_error(
    placed_coords,
    crosslinker,
    candidate_group,
    path_coords,
    crosslink_bond_length,
    pbc,
    box_lengths,
):
    """Compute the maximum fractional bond length error for a placement.

    Checks each connection site against each backbone bead in its group.

    Returns
    -------
    max_error : float
        Maximum |actual_dist - target| / target across all bonds.
    """
    conn_sites = crosslinker.connection_sites
    max_error = 0.0

    for conn_idx, phys_node in enumerate(conn_sites):
        site_pos = placed_coords[phys_node]
        for bb_idx in candidate_group[conn_idx]:
            dist = float(_pbc_distance(site_pos, path_coords[bb_idx], box_lengths, pbc))
            error = abs(dist - crosslink_bond_length) / crosslink_bond_length
            max_error = max(max_error, error)

    return max_error


# =============================================================================
# Placement with Rotation Search
# =============================================================================


def _try_placement_with_rotations(
    crosslinker,
    candidate_group,
    path_coords,
    crosslink_bond_length,
    excluded_indices,
    pbc,
    box_lengths,
    tolerance,
    minimum_separation,
    n_rotation_samples,
    rng,
):
    """Attempt to place crosslinker with rotation search for overlap avoidance.

    1. Computes the optimal placement (translation + alignment).
    2. If single-site or the base placement works, returns it.
    3. Otherwise, tries rotations around the centroid to avoid overlaps.

    Returns
    -------
    placed_coords : np.ndarray or None
        Final placed coordinates, or None if no valid placement found.
    """
    # Compute base placement
    base_coords, target_node_positions = _compute_optimal_placement(
        crosslinker,
        candidate_group,
        path_coords,
        crosslink_bond_length,
        pbc,
        box_lengths,
        tolerance,
    )

    # Check bond lengths on base placement
    base_error = _compute_max_bond_error(
        base_coords,
        crosslinker,
        candidate_group,
        path_coords,
        crosslink_bond_length,
        pbc,
        box_lengths,
    )

    if base_error > tolerance:
        # Base placement doesn't satisfy bond lengths — skip
        return None

    # Check overlaps on base placement
    if not _has_overlaps(
        base_coords, path_coords, excluded_indices, minimum_separation, pbc, box_lengths
    ):
        return base_coords

    # --- Need rotational search to avoid overlaps ---
    # For single-site crosslinkers, rotation doesn't help
    if crosslinker.n_sites == 1:
        return None

    # Determine rotation axis
    unique_nodes = list(target_node_positions.keys())
    centroid = base_coords.mean(axis=0)

    if len(unique_nodes) >= 2:
        # Rotate around axis connecting the physical nodes
        positions = np.array([target_node_positions[n] for n in unique_nodes])
        axis = positions[-1] - positions[0]
    elif len(unique_nodes) == 1:
        # Rotate around the axis from centroid to the single physical node
        axis = target_node_positions[unique_nodes[0]] - centroid
    else:
        axis = np.array([0, 0, 1], dtype=np.float32)

    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-10:
        axis = np.array([0, 0, 1], dtype=np.float32)
    else:
        axis = axis / axis_norm

    # Rotation pivot: the physical node(s) must stay fixed
    # Rotate around the axis through the centroid of physical nodes
    if len(unique_nodes) >= 1:
        pivot = np.mean([target_node_positions[n] for n in unique_nodes], axis=0)
    else:
        pivot = centroid

    best_coords = None
    best_error = float("inf")

    angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)
    rng.shuffle(angles)

    for angle in angles:
        # Rotate around axis through pivot
        relative = base_coords - pivot
        rotated = _rotate_around_axis(relative, axis, angle) + pivot

        # Verify bond lengths still valid after rotation
        error = _compute_max_bond_error(
            rotated,
            crosslinker,
            candidate_group,
            path_coords,
            crosslink_bond_length,
            pbc,
            box_lengths,
        )
        if error > tolerance:
            continue

        # Check overlaps
        if not _has_overlaps(
            rotated, path_coords, excluded_indices, minimum_separation, pbc, box_lengths
        ):
            if error < best_error:
                best_error = error
                best_coords = rotated.copy()
                break  # First overlap-free is good enough

    return best_coords


# =============================================================================
# Path Insertion
# =============================================================================


def _insert_crosslinker(path, crosslinker, placed_coords, candidate_group):
    """Insert crosslinker into path: add nodes, internal edges, and external bonds.

    Parameters
    ----------
    path : Path
    crosslinker : CrosslinkerGeometry
    placed_coords : np.ndarray, shape (n_sites, 3)
    candidate_group : list of list[int]
        candidate_group[conn_idx] = list of backbone node indices for that
        connection site.
    """
    n_existing = len(path.coordinates)
    n_new = crosslinker.n_sites
    conn_sites = crosslinker.connection_sites

    # Append coordinates and bead names
    path.coordinates = np.vstack([path.coordinates, placed_coords])
    path.beads = np.append(path.beads, crosslinker.beads)

    # Add nodes
    new_indices = list(range(n_existing, n_existing + n_new))
    for i, new_idx in enumerate(new_indices):
        path.bond_graph.add_node(
            new_idx,
            name=str(crosslinker.beads[i]),
            xyz=placed_coords[i],
        )

    # Add internal crosslinker edges
    for u, v in crosslinker.bond_graph.edges():
        bond_type = (str(crosslinker.beads[u]), str(crosslinker.beads[v]))
        path.bond_graph.add_edge(new_indices[u], new_indices[v], bond_type=bond_type)

    # Add external bonds: each connection site bonds to its backbone group
    for conn_idx, phys_node in enumerate(conn_sites):
        new_conn_node = new_indices[phys_node]
        for bb_node in candidate_group[conn_idx]:
            # Avoid duplicate edges (same physical node may appear multiple times)
            if not path.bond_graph.has_edge(new_conn_node, bb_node):
                bond_type = (
                    str(crosslinker.beads[phys_node]),
                    str(path.beads[bb_node]),
                )
                path.bond_graph.add_edge(new_conn_node, bb_node, bond_type=bond_type)


# =============================================================================
# Main Crosslink Function
# =============================================================================


def crosslink(
    path,
    crosslinker=None,
    backbone_name="_A",
    crosslink_bond_length=0.2,
    max_backbone_degree=4,
    tolerance=0.1,
    excluded_bond_depth=2,
    volume_constraint=None,
    seed=42,
    n_rotation_samples=36,
    minimum_separation=None,
):
    """Place a crosslinker bonded to backbone beads, modifying path in place.

    Parameters
    ----------
    path : Path
        The Path object to modify in place.
    crosslinker : CrosslinkerGeometry, optional
        If None, auto-generated from bead_name and n_connection_sites.
    backbone_name : str or tuple
        Specifies what each connection site bonds to:
        - String: every connection site bonds to one bead of that type.
        - Tuple of length == n_connection_sites:
            - String entry: bonds to one bead.
            - Tuple entry (e.g., ("_B", "_B")): bonds to multiple neighboring beads.

        Examples::

            backbone_name = "_B"
            # All connection sites bond to one "_B" each.

            backbone_name = (("_B", "_B"), "_B")
            # Site 0 bonds to two neighboring "_B" beads.
            # Site 1 bonds to one "_B" bead.

            backbone_name = (("_B", "_B"), ("_B", "_B"))
            # Both sites bond to two neighboring "_B" beads each.
            # If connection_sites=[0,0], ALL four beads must be at
            # crosslink_bond_length from the same physical node.

    crosslink_bond_length : float, default 0.2
        Desired distance from crosslinker connection node to each backbone bead.
    max_backbone_degree : int, default 4
        Maximum graph degree a backbone node may have to remain eligible.
    tolerance : float, default 0.1
        Fractional tolerance on bond lengths (0.1 = ±10%).
    excluded_bond_depth : int, default 2
        Beads within this many bonds of selected nodes are excluded from
        being selected as additional connection points.
    volume_constraint : optional
        Provides PBC information.
    seed : int, default 42
    n_rotation_samples : int, default 36
    minimum_separation : float, optional
        Minimum distance to existing beads.
        Defaults to crosslink_bond_length / 4.
        This will only apply to non-bonded beads,
        i.e. you can set a separation larger than
        1/2 the crosslink_bond_length and still
        achieve valid placements.

    Returns
    -------
    path : Path (modified in place)
        Path object with a single added crosslink.
    """
    rng = np.random.default_rng(seed + len(path.coordinates))

    # --- Default crosslinker ---
    if crosslinker is None:  # Single two site crosslinker
        crosslinker = CrosslinkerGeometry.single_site()

    # --- Parse backbone specification ---
    backbone_specs = _parse_backbone_spec(backbone_name, crosslinker.connection_sites)

    # --- PBC setup ---
    pbc, box_lengths = _get_pbc_info(volume_constraint)

    # --- Min separation ---
    if minimum_separation is None:
        minimum_separation = crosslink_bond_length * 0.25

    # --- Find candidate backbone groups ---
    candidate_groups = _find_candidate_groups(
        path,
        crosslinker,
        backbone_specs,
        crosslink_bond_length,
        tolerance,
        excluded_bond_depth,
        max_backbone_degree,
        pbc,
        box_lengths,
        rng,
    )

    if not candidate_groups:
        n_beads_in_cl = len(
            crosslinker.coordinates
        )  # divisor to approx previous number of crosslinks
        n_clinks = (
            sum(1 for b in path.beads if b in set(crosslinker.beads)) // n_beads_in_cl
        )
        raise PathConvergenceError(
            f"Could not find backbone beads matching crosslinker geometry with "
            f"bond_length={crosslink_bond_length} (tolerance=±{tolerance * 100:.0f}%).\n"
            f"max_backbone_degree={max_backbone_degree}, "
            f"excluded_bond_depth={excluded_bond_depth}.\n"
            f"Current crosslinks: {n_clinks}.\n"
            "Ways to increase crosslinking:\n"
            " - Increase crosslink_bond_length or tolerance\n"
            " - Increase max_backbone_degree\n"
            " - Decrease excluded_bond_depth\n"
            " - Pack at higher density\n"
        )

    # --- Try each candidate group ---
    rng.shuffle(candidate_groups)
    placed = False

    for candidate_group in candidate_groups:
        # Collect all backbone nodes involved
        all_backbone_nodes = _all_beads_from_groups(candidate_group)

        # Excluded indices for overlap check: backbone nodes + their immediate neighbors
        excluded_indices = set(all_backbone_nodes)
        for node in all_backbone_nodes:
            if node in path.bond_graph:
                for neighbor in path.bond_graph.neighbors(node):
                    excluded_indices.add(neighbor)

        # Attempt placement
        placed_coords = _try_placement_with_rotations(
            crosslinker,
            candidate_group,
            path.coordinates,
            crosslink_bond_length,
            excluded_indices,
            pbc,
            box_lengths,
            tolerance,
            minimum_separation,
            n_rotation_samples,
            rng,
        )

        if placed_coords is None:
            continue

        # --- Success: insert into path ---
        _insert_crosslinker(path, crosslinker, placed_coords, candidate_group)

        if volume_constraint is not None:
            path.wrap_inside_box(volume_constraint)

        placed = True
        break

    if not placed:
        n_clinks = sum(1 for b in path.beads if b in set(crosslinker.beads))
        raise PathConvergenceError(
            f"Could not place crosslinker without overlaps. "
            f"Tried {len(candidate_groups)} candidate groups.\n"
            f"Current crosslinks: {n_clinks}.\n"
            "Consider increasing tolerance, n_rotation_samples, or decreasing minimum_separation."
        )

    return path


# =============================================================================
# PBC Info Extraction
# =============================================================================


def _get_pbc_info(volume_constraint):
    """Extract PBC flags and box lengths from a volume constraint.

    Returns
    -------
    pbc : np.ndarray, shape (3,), dtype bool
    box_lengths : np.ndarray, shape (3,), dtype float32
    """
    if volume_constraint is None:
        return np.array([False, False, False], dtype=bool), np.zeros(
            3, dtype=np.float32
        )

    if isinstance(volume_constraint, CuboidConstraint):
        pbc = np.array(volume_constraint.pbc, dtype=bool)
        box_lengths = np.array(volume_constraint.box_lengths, dtype=np.float32)
    elif isinstance(volume_constraint, CylinderConstraint):
        pbc = np.array(
            [False, False, getattr(volume_constraint, "pbc_z", False)], dtype=bool
        )
        box_lengths = np.array(
            [0.0, 0.0, getattr(volume_constraint, "length", 0.0)], dtype=np.float32
        )
    elif hasattr(volume_constraint, "pbc"):
        pbc = np.array(volume_constraint.pbc, dtype=bool)
        box_lengths = np.array(
            [
                getattr(volume_constraint, "Lx", 0),
                getattr(volume_constraint, "Ly", 0),
                getattr(volume_constraint, "Lz", 0),
            ],
            dtype=np.float32,
        )
    else:
        pbc = np.array([False, False, False], dtype=bool)
        box_lengths = np.zeros(3, dtype=np.float32)

    return pbc, box_lengths
