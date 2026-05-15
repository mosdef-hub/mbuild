from itertools import permutations

import networkx as nx
import numpy as np

from mbuild.exceptions import PathConvergenceError
from mbuild.path.build import Path
from mbuild.path.constraints import CuboidConstraint, CylinderConstraint
from mbuild.path.path_utils import calculate_sq_distances


class CrosslinkerGeometry(Path):
    """A Path subclass representing a crosslinker with connection site metadata.

    Inherits all Path functionality (coordinates, bond_graph, beads) and adds
    metadata about which sites form bonds to external beads.

    Parameters
    ----------
    coordinates : array-like, shape (N, 3)
        Positions of crosslinker sites. Will be re-centered at origin.
    bond_graph : networkx.Graph, optional
        Internal connectivity.
    bead_name : str or np.ndarray, default "_R"
        Name(s) for each site.
    connection_sites : list of int
        Node indices that form bonds to external beads.
        Length determines how many external bonds this geometry requires.
        Repeated indices allow one site to bond to multiple external beads.

    Examples
    --------
    >>> cl = CrosslinkerGeometry.equilateral_triangle(bond_length=0.27)
    >>> cl.connection_sites
    [0, 1, 2]
    >>> cl.bond_graph.edges()
    EdgeView([(0, 1), (0, 2), (1, 2)])
    """

    def __init__(
        self,
        coordinates=None,
        bond_graph=None,
        bead_name="_R",
        connection_sites=None,
    ):
        super().__init__(
            coordinates=coordinates, bond_graph=bond_graph, bead_name=bead_name
        )

        if connection_sites is None:
            connection_sites = list(range(len(self.coordinates)))

        self.connection_sites = list(connection_sites)
        self._validate_connection_sites()

        # Center positions at origin
        if len(self.coordinates) > 0:
            centroid = np.mean(self.coordinates, axis=0)
            self.coordinates = (self.coordinates - centroid).astype(np.float32)

    def _validate_connection_sites(self):
        n = len(self.coordinates)
        for idx in self.connection_sites:
            if not (0 <= idx < n):
                raise ValueError(
                    f"Connection site {idx} is out of range (must be 0 to {n - 1})"
                )

    @property
    def n_sites(self) -> int:
        return len(self.coordinates)

    @property
    def n_connections(self) -> int:
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
    def equilateral_triangle(
        cls, bond_length=0.27, bead_name="_R", connection_sites=None
    ):
        """Three sites in a flat equilateral triangle."""
        R = bond_length / np.sqrt(3)
        angles = [0, 2 * np.pi / 3, 4 * np.pi / 3]
        positions = np.array(
            [[R * np.cos(a), R * np.sin(a), 0.0] for a in angles], dtype=np.float32
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
# Orientation / Rotation Utilities
# =============================================================================


def _kabsch_rotation(P, Q):
    """Optimal rotation aligning P onto Q (Kabsch algorithm)."""
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, d])
    R = Vt.T @ sign_matrix @ U.T
    return R.astype(np.float32)


def _rotation_between_vectors(a, b):
    """Rotation matrix rotating unit vector a onto b (Rodrigues)."""
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = float(np.linalg.norm(v))

    if s < 1e-10:
        if c > 0:
            return np.eye(3, dtype=np.float32)
        perp = np.array([1, 0, 0], dtype=np.float32)
        if abs(np.dot(a, perp)) > 0.9:
            perp = np.array([0, 1, 0], dtype=np.float32)
        perp = perp - np.dot(perp, a) * a
        perp /= np.linalg.norm(perp)
        return (2 * np.outer(perp, perp) - np.eye(3)).astype(np.float32)

    vx = np.array(
        [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype=np.float32
    )
    R = np.eye(3, dtype=np.float32) + vx + vx @ vx * ((1 - c) / (s * s))
    return R


def _random_rotation_matrix(rng):
    """Uniform random rotation via QR decomposition."""
    H = rng.standard_normal((3, 3)).astype(np.float32)
    Q, R = np.linalg.qr(H)
    Q *= np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q.astype(np.float32)


def _rotation_about_axis(axis, angle):
    """Rotation matrix for rotation by `angle` radians about `axis`."""
    axis = np.asarray(axis, dtype=np.float32)
    axis = axis / np.linalg.norm(axis)
    K = np.array(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ],
        dtype=np.float32,
    )
    R = np.eye(3, dtype=np.float32) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    return R


def _orient_replacement(
    replacement,
    neighbor_vectors,
    nearby_coords=None,
    center=None,
    rng=None,
    n_rotation_samples=36,
):
    """Orient a replacement geometry so connection sites point toward neighbors.

    Performs Kabsch alignment then searches residual rotational DOF to
    minimize overlaps with nearby non-bonded sites.

    Parameters
    ----------
    replacement : CrosslinkerGeometry
        The replacement structure (centered at origin).
    neighbor_vectors : np.ndarray, shape (n_connections, 3)
        Vectors from placement center to each neighbor being bonded.
    nearby_coords : np.ndarray, shape (M, 3), optional
        Coordinates of nearby non-bonded sites (for overlap minimization).
        Should be relative to the placement center.
    center : np.ndarray, shape (3,), optional
        The absolute position where replacement is centered (for overlap calc).
    rng : numpy.random.Generator, optional
    n_rotation_samples : int, default 36
        Number of rotational samples for overlap minimization.

    Returns
    -------
    oriented_positions : np.ndarray, shape (n_sites, 3)
        Rotated replacement positions (still centered at origin).
    assignment : list of int
        assignment[i] = j means connection_sites[i] bonds to neighbor j.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    positions = replacement.coordinates.copy()
    unique_conn = replacement.unique_connection_sites
    conn_positions = positions[unique_conn]

    n_unique = len(unique_conn)

    # --- Degenerate case: all connection sites at origin ---
    if np.allclose(conn_positions, 0, atol=1e-8):
        R = _random_rotation_matrix(rng)
        oriented = (R @ positions.T).T
        return oriented, list(range(replacement.n_connections))

    # --- Group target vectors by unique connection site (for repeated sites) ---
    conn_site_to_targets = {}
    for i, site_idx in enumerate(replacement.connection_sites):
        conn_site_to_targets.setdefault(site_idx, []).append(i)

    # Use centroid of targets for each unique connection site
    unique_targets = np.array(
        [
            np.mean(neighbor_vectors[conn_site_to_targets[s]], axis=0)
            for s in unique_conn
        ],
        dtype=np.float32,
    )

    # --- Normalize for angular alignment ---
    P_norms = np.linalg.norm(conn_positions, axis=1, keepdims=True)
    P_normed = conn_positions / np.maximum(P_norms, 1e-10)

    T_norms = np.linalg.norm(unique_targets, axis=1, keepdims=True)
    T_normed = unique_targets / np.maximum(T_norms, 1e-10)

    # --- Single connection site case ---
    if n_unique == 1:
        R = _rotation_between_vectors(P_normed[0], T_normed[0])
        best_R = R
        best_perm = [0]
    else:
        # --- Try all permutations (feasible for n <= 5) ---
        best_perm = list(range(n_unique))
        best_R = np.eye(3, dtype=np.float32)
        best_cost = np.inf

        for perm in permutations(range(n_unique)):
            Q_perm = T_normed[list(perm)]
            R = _kabsch_rotation(P_normed, Q_perm)
            rotated = (R @ P_normed.T).T
            cost = np.sum((rotated - Q_perm) ** 2)
            if cost < best_cost:
                best_cost = cost
                best_perm = list(perm)
                best_R = R

    oriented = (best_R @ positions.T).T

    # --- Overlap minimization: search residual rotation ---
    if nearby_coords is not None and len(nearby_coords) > 0 and n_rotation_samples > 1:
        # Determine a rotation axis (mean direction of connection sites after alignment)
        rotated_conn = oriented[unique_conn]
        mean_axis = np.mean(rotated_conn, axis=0)
        mean_axis_norm = np.linalg.norm(mean_axis)

        if mean_axis_norm > 1e-8:
            rotation_axis = mean_axis / mean_axis_norm
        else:
            # Use normal to the plane of connection sites or random axis
            rotation_axis = np.array([0, 0, 1], dtype=np.float32)

        best_oriented = oriented.copy()
        best_min_dist = _min_distance_to_neighbors(oriented, nearby_coords)

        angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)
        for angle in angles:
            R_trial = _rotation_about_axis(rotation_axis, angle)
            trial = (R_trial @ oriented.T).T
            min_dist = _min_distance_to_neighbors(trial, nearby_coords)
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_oriented = trial

        oriented = best_oriented

    # --- Build full assignment ---
    full_assignment = _build_full_assignment(
        replacement, best_perm, unique_conn, conn_site_to_targets
    )

    return oriented, full_assignment


def _orient_replacement_partial(
    replacement,
    partial_conn_sites,
    neighbor_vectors,
    nearby_coords=None,
    center=None,
    rng=None,
    n_rotation_samples=36,
):
    """Orient a replacement when fewer neighbors exist than connection sites.

    Uses only the first ``len(neighbor_vectors)`` connection sites for alignment.
    The remaining connection sites are left dangling (not bonded).

    Parameters
    ----------
    replacement : CrosslinkerGeometry
        The replacement structure (centered at origin).
    partial_conn_sites : list of int
        Subset of replacement.connection_sites to use for alignment.
        Length matches len(neighbor_vectors).
    neighbor_vectors : np.ndarray, shape (degree, 3)
        Vectors from placement center to each neighbor being bonded.
    nearby_coords : np.ndarray, shape (M, 3), optional
        Coordinates of nearby non-bonded sites for overlap minimization.
    center : np.ndarray, shape (3,), optional
        Absolute position of placement center.
    rng : numpy.random.Generator, optional
    n_rotation_samples : int, default 36
        Number of rotational samples for overlap minimization.

    Returns
    -------
    oriented_positions : np.ndarray, shape (n_sites, 3)
        Rotated replacement positions (still centered at origin).
    assignment : list of int
        assignment[i] = j means partial_conn_sites[i] bonds to neighbor j.
        Length equals len(partial_conn_sites) (i.e., degree).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    positions = replacement.coordinates.copy()
    degree = len(neighbor_vectors)

    if degree == 0:
        R = _random_rotation_matrix(rng)
        return (R @ positions.T).T, []

    # Get positions of the partial connection sites
    unique_partial = list(dict.fromkeys(partial_conn_sites))
    conn_positions = positions[unique_partial]

    # If all connection sites at origin, random rotation
    if np.allclose(conn_positions, 0, atol=1e-8):
        R = _random_rotation_matrix(rng)
        oriented = (R @ positions.T).T
        return oriented, list(range(degree))

    # Group partial connections by unique site
    conn_site_to_targets = {}
    for i, site_idx in enumerate(partial_conn_sites):
        conn_site_to_targets.setdefault(site_idx, []).append(i)

    # Use centroid of targets for each unique connection site
    unique_targets = np.array(
        [
            np.mean(neighbor_vectors[conn_site_to_targets[s]], axis=0)
            for s in unique_partial
        ],
        dtype=np.float32,
    )

    # Normalize for angular alignment
    P_norms = np.linalg.norm(conn_positions, axis=1, keepdims=True)
    P_normed = conn_positions / np.maximum(P_norms, 1e-10)

    T_norms = np.linalg.norm(unique_targets, axis=1, keepdims=True)
    T_normed = unique_targets / np.maximum(T_norms, 1e-10)

    n_unique = len(unique_partial)

    if n_unique == 1:
        R = _rotation_between_vectors(P_normed[0], T_normed[0])
        best_R = R
        best_perm = [0]
    else:
        # Try all permutations of unique sites -> unique targets
        best_perm = list(range(n_unique))
        best_R = np.eye(3, dtype=np.float32)
        best_cost = np.inf

        for perm in permutations(range(n_unique)):
            Q_perm = T_normed[list(perm)]
            R = _kabsch_rotation(P_normed, Q_perm)
            rotated = (R @ P_normed.T).T
            cost = np.sum((rotated - Q_perm) ** 2)
            if cost < best_cost:
                best_cost = cost
                best_perm = list(perm)
                best_R = R

    oriented = (best_R @ positions.T).T

    # --- Overlap minimization ---
    if nearby_coords is not None and len(nearby_coords) > 0 and n_rotation_samples > 1:
        rotated_conn = oriented[unique_partial]
        mean_axis = np.mean(rotated_conn, axis=0)
        mean_axis_norm = np.linalg.norm(mean_axis)

        if mean_axis_norm > 1e-8:
            rotation_axis = mean_axis / mean_axis_norm
        else:
            rotation_axis = np.array([0, 0, 1], dtype=np.float32)

        best_oriented = oriented.copy()
        best_min_dist = _min_distance_to_neighbors(oriented, nearby_coords)

        angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)
        for angle in angles:
            R_trial = _rotation_about_axis(rotation_axis, angle)
            trial = (R_trial @ oriented.T).T
            min_dist = _min_distance_to_neighbors(trial, nearby_coords)
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_oriented = trial

        oriented = best_oriented

    # Build assignment: maps conn_idx (0..degree-1) -> neighbor_idx
    full_assignment = _build_full_assignment_partial(
        partial_conn_sites, best_perm, unique_partial, conn_site_to_targets
    )

    return oriented, full_assignment


def _build_full_assignment_partial(
    partial_conn_sites, perm, unique_partial, conn_site_to_targets
):
    """Build assignment for partial connection case.

    Parameters
    ----------
    partial_conn_sites : list of int
        The subset of connection sites being used.
    perm : list of int
        Permutation of unique sites to unique targets.
    unique_partial : list of int
        Unique site indices in partial_conn_sites.
    conn_site_to_targets : dict
        Maps site_idx -> list of connection indices using that site.

    Returns
    -------
    full_assignment : list of int
        assignment[i] = neighbor_index for partial_conn_sites[i].
    """
    degree = len(partial_conn_sites)

    if degree == len(unique_partial) and all(
        len(v) == 1 for v in conn_site_to_targets.values()
    ):
        # Simple case: no repeated sites, 1:1 mapping
        return perm

    # General case with possible repeated sites
    full_assignment = [0] * degree
    site_counter = {}

    for conn_idx, site_idx in enumerate(partial_conn_sites):
        unique_pos = unique_partial.index(site_idx)
        target_neighbor_group = perm[unique_pos]

        count = site_counter.get(site_idx, 0)
        target_site = unique_partial[target_neighbor_group]
        target_indices = conn_site_to_targets[target_site]
        full_assignment[conn_idx] = target_indices[count % len(target_indices)]
        site_counter[site_idx] = count + 1

    return full_assignment


def _min_distance_to_neighbors(positions, nearby_coords):
    """Compute minimum distance between any replacement site and nearby sites."""
    if len(nearby_coords) == 0:
        return np.inf
    # Pairwise distances
    diffs = positions[:, None, :] - nearby_coords[None, :, :]
    dists = np.linalg.norm(diffs, axis=-1)
    return np.min(dists)


def _build_full_assignment(replacement, perm, unique_conn, conn_site_to_targets):
    """Map each connection_sites entry to a neighbor index using the permutation.

    Parameters
    ----------
    replacement : CrosslinkerGeometry
    perm : list of int
        Permutation mapping unique connection site positions to neighbor group positions.
    unique_conn : list of int
        Unique connection site indices.
    conn_site_to_targets : dict
        Maps site_idx -> list of connection indices using that site.

    Returns
    -------
    full_assignment : list of int
        assignment[i] = neighbor_index for connection_sites[i].
    """
    n_connections = replacement.n_connections

    if n_connections == len(unique_conn):
        # No repeated sites: perm directly maps connection -> neighbor
        return perm

    # With repeated sites: need to assign individual connections within groups
    full_assignment = [0] * n_connections

    # perm[unique_pos] tells us which neighbor group this unique site maps to
    # We need to iterate through connection_sites and assign appropriately
    site_counter = {}  # track how many times each site has been assigned
    for conn_idx, site_idx in enumerate(replacement.connection_sites):
        unique_pos = unique_conn.index(site_idx)
        target_neighbor_group = perm[unique_pos]

        # Round-robin within the target group
        count = site_counter.get(site_idx, 0)
        # The target_neighbor_group-th unique site's target list
        target_site = unique_conn[target_neighbor_group]
        target_indices = conn_site_to_targets[target_site]
        full_assignment[conn_idx] = target_indices[count % len(target_indices)]
        site_counter[site_idx] = count + 1

    return full_assignment


# =============================================================================
# General replace_sites function
# =============================================================================


def replace_sites(
    path,
    replacement,
    sites=None,
    bead_name=None,
    bond_length=None,
    tolerance=0.1,
    min_separation=None,
    volume_constraint=None,
    n_rotation_samples=36,
    overlap_radius=None,
    seed=42,
):
    """Replace sites in a Path with a multi-site substructure, in place.

    Uses rigid-body alignment to place the replacement such that:
    - Connection sites are at specified bond_length from their neighbors
    - Internal geometry is preserved exactly (rigid body)
    - Placed sites are at least min_separation from non-bonded neighbors

    The optimization sweeps residual rotational degrees of freedom and, if
    needed, offsets the placement perpendicular to the constraint axis to
    find a clash-free configuration.

    Parameters
    ----------
    path : Path
        The path to modify in place.
    replacement : CrosslinkerGeometry
        The replacement structure.
    sites : list of int, optional
        Explicit node indices to replace. Mutually exclusive with ``bead_name``.
    bead_name : str, optional
        Replace all sites matching this bead name. Mutually exclusive with ``sites``.
    bond_length : float, optional
        Desired bond length from each connection site to its neighbor.
        If None, uses the existing distance from the replaced site to each neighbor.
    tolerance : float, default 0.1
        Fractional tolerance on bond lengths (±10%).
    min_separation : float, optional
        Minimum allowed distance between any placed crosslinker site and any
        existing non-bonded site. If None, no overlap checking is performed.
    volume_constraint : CuboidConstraint or CylinderConstraint, optional
        For periodic boundary handling.
    n_rotation_samples : int, default 36
        Angular samples for residual DOF optimization.
    overlap_radius : float, optional
        Radius around each replaced site to gather nearby sites for overlap check.
        If None, defaults to 3x max extent of replacement geometry.
    seed : int, default 42
        Random seed.

    Returns
    -------
    Path
        The same path object (modified in place), for method chaining.
    """
    if sites is None and bead_name is None:
        raise ValueError("Must specify either `sites` or `bead_name`.")
    if sites is not None and bead_name is not None:
        raise ValueError("Cannot specify both `sites` and `bead_name`.")

    rng = np.random.default_rng(seed)

    if bead_name is not None:
        sites = [
            node for node in path.bond_graph.nodes() if path.beads[node] == bead_name
        ]
    sites = list(sites)

    if len(sites) == 0:
        return path

    for site in sites:
        degree = path.bond_graph.degree[site]
        if degree > replacement.n_connections:
            raise ValueError(
                f"Site {site} has degree {degree} but replacement only has "
                f"{replacement.n_connections} connection sites. "
                f"Degree must not exceed len(replacement.connection_sites)."
            )

    pbc, box_lengths = _get_pbc_info(volume_constraint)

    if overlap_radius is None:
        if replacement.n_sites > 0:
            max_extent = np.max(np.linalg.norm(replacement.coordinates, axis=1))
            overlap_radius = max(max_extent * 3, 0.5)
        else:
            overlap_radius = 0.5

    sites_set = set(sites)

    # --- Compute all replacements ---
    replacement_data = {}

    for site in sites:
        neighbors = list(path.bond_graph.neighbors(site))
        site_pos = path.coordinates[site]
        degree = len(neighbors)

        # Unwrapped neighbor positions (absolute)
        neighbor_coords = np.zeros((degree, 3), dtype=np.float32)
        for i, nb in enumerate(neighbors):
            delta = path.coordinates[nb] - site_pos
            delta = _apply_mic(delta, pbc, box_lengths)
            neighbor_coords[i] = site_pos + delta

        # Target bond lengths
        if bond_length is not None:
            target_bond_lengths = np.full(degree, bond_length, dtype=np.float32)
        else:
            target_bond_lengths = np.array(
                [np.linalg.norm(neighbor_coords[i] - site_pos) for i in range(degree)],
                dtype=np.float32,
            )

        # Nearby non-bonded coords (absolute positions for overlap checking)
        nearby_coords = _get_nearby_nonbonded_absolute(
            path, site, neighbors, site_pos, overlap_radius, pbc, box_lengths, sites_set
        )

        # Active connection sites (may be fewer than total if degree < n_connections)
        active_conn_sites = replacement.connection_sites[:degree]

        # Place via rigid body alignment + overlap optimization
        placed_positions, assignment = _place_replacement_rigid(
            replacement,
            active_conn_sites,
            neighbor_coords,
            target_bond_lengths,
            nearby_coords,
            rng,
            n_rotation_samples,
            tolerance,
            min_separation,
        )

        # Wrap into PBC
        placed_positions = _wrap_positions(placed_positions, pbc, box_lengths)

        replacement_data[site] = {
            "positions": placed_positions,
            "neighbors": neighbors,
            "assignment": assignment,
        }

    # --- Rebuild path in place ---
    kept_nodes = [n for n in sorted(path.bond_graph.nodes()) if n not in sites_set]

    new_coords = []
    new_beads = []
    old_to_new = {}
    current_idx = 0

    for old_node in kept_nodes:
        old_to_new[old_node] = current_idx
        new_coords.append(path.coordinates[old_node])
        new_beads.append(path.beads[old_node])
        current_idx += 1

    site_new_indices = {}
    for site in sites:
        data = replacement_data[site]
        indices = []
        for i in range(replacement.n_sites):
            indices.append(current_idx)
            new_coords.append(data["positions"][i])
            new_beads.append(replacement.beads[i])
            current_idx += 1
        site_new_indices[site] = indices

    # --- Build edges ---
    new_edges = []

    # Kept-to-kept
    for u, v, data in path.bond_graph.edges(data=True):
        if u in sites_set or v in sites_set:
            continue
        new_edges.append((old_to_new[u], old_to_new[v], data))

    # Internal replacement bonds
    for site in sites:
        indices = site_new_indices[site]
        for i, j in replacement.bond_graph.edges():
            edge_data = {
                "bond_type": (str(replacement.beads[i]), str(replacement.beads[j])),
            }
            new_edges.append((indices[i], indices[j], edge_data))

    # Connection-to-neighbor bonds
    for site in sites:
        data = replacement_data[site]
        indices = site_new_indices[site]
        neighbors = data["neighbors"]
        assignment = data["assignment"]

        for conn_idx, neighbor_idx in enumerate(assignment):
            repl_site_local = replacement.connection_sites[conn_idx]
            repl_new_node = indices[repl_site_local]

            external_neighbor = neighbors[neighbor_idx]

            if external_neighbor in old_to_new:
                ext_new_node = old_to_new[external_neighbor]
                ext_bead_name = path.beads[external_neighbor]
            elif external_neighbor in site_new_indices:
                other_data = replacement_data[external_neighbor]
                other_neighbors = other_data["neighbors"]
                try:
                    back_idx = other_neighbors.index(site)
                except ValueError:
                    continue
                other_assignment = other_data["assignment"]
                found = False
                for other_conn_idx, other_target in enumerate(other_assignment):
                    if other_target == back_idx:
                        other_repl_site_local = replacement.connection_sites[
                            other_conn_idx
                        ]
                        ext_new_node = site_new_indices[external_neighbor][
                            other_repl_site_local
                        ]
                        ext_bead_name = replacement.beads[other_repl_site_local]
                        found = True
                        break
                if not found:
                    continue
            else:
                continue

            edge_data = {
                "bond_type": (
                    str(replacement.beads[repl_site_local]),
                    str(ext_bead_name),
                ),
            }
            new_edges.append((repl_new_node, ext_new_node, edge_data))

    # --- Overwrite path ---
    path.coordinates = np.array(new_coords, dtype=np.float32)
    path.beads = np.array(new_beads, dtype="U10")

    path.bond_graph = nx.Graph()
    for i in range(len(path.coordinates)):
        path.bond_graph.add_node(i, name=str(path.beads[i]), xyz=path.coordinates[i])

    seen_edges = set()
    for u, v, data in new_edges:
        edge_key = (min(u, v), max(u, v))
        if edge_key not in seen_edges:
            seen_edges.add(edge_key)
            path.bond_graph.add_edge(u, v, **data)

    return path


def _get_nearby_nonbonded_absolute(
    path, site, neighbors, site_pos, radius, pbc, box_lengths, excluded_sites
):
    """Get absolute positions of nearby non-bonded sites for overlap checking.

    Unlike the relative version, returns absolute (unwrapped relative to site_pos)
    positions suitable for direct distance comparison with placed crosslinker sites.

    Parameters
    ----------
    path : Path
    site : int
    neighbors : list of int
    site_pos : np.ndarray, shape (3,)
    radius : float
    pbc : array-like of bool
    box_lengths : np.ndarray
    excluded_sites : set of int

    Returns
    -------
    np.ndarray, shape (M, 3)
    """
    neighbors_set = set(neighbors) | {site} | excluded_sites

    nearby = []
    for node in path.bond_graph.nodes():
        if node in neighbors_set:
            continue
        delta = path.coordinates[node] - site_pos
        delta = _apply_mic(delta, pbc, box_lengths)
        dist = np.linalg.norm(delta)
        if dist <= radius:
            # Return absolute unwrapped position
            nearby.append(site_pos + delta)

    if len(nearby) == 0:
        return np.empty((0, 3), dtype=np.float32)
    return np.array(nearby, dtype=np.float32)


def _place_replacement_rigid(
    replacement,
    active_conn_sites,
    neighbor_coords,
    target_bond_lengths,
    nearby_coords,
    rng,
    n_rotation_samples,
    tolerance,
    min_separation=None,
):
    """Place replacement via rigid-body alignment guaranteeing bond lengths and no overlaps.

    Strategy:
    1. Compute target positions for connection sites (at bond_length from neighbors).
    2. Solve Procrustes alignment (rotation + translation) to map connection sites to targets.
    3. For residual rotational DOF (1 DOF for 2 unique sites, full DOF for 1 unique site):
       - Sweep angles and pick the one that:
         a) Keeps all connection-to-neighbor distances within tolerance
         b) Maximizes minimum distance to nearby non-bonded sites (avoids overlaps)
    4. If no rotation satisfies min_separation, offset the center along the rotation axis
       (push the crosslinker "above" or "below" the backbone plane) and re-sweep.
    5. As a last resort, return the best-scoring rotation even if below min_separation.

    Parameters
    ----------
    replacement : CrosslinkerGeometry
        The replacement structure (centered at origin).
    active_conn_sites : list of int
        Subset of replacement.connection_sites to bond (length = degree).
    neighbor_coords : np.ndarray, shape (degree, 3)
        Unwrapped absolute positions of neighboring beads.
    target_bond_lengths : np.ndarray, shape (degree,)
        Desired bond length from each connection site to its neighbor.
    nearby_coords : np.ndarray, shape (M, 3)
        Absolute positions of nearby non-bonded sites for overlap avoidance.
    rng : numpy.random.Generator
    n_rotation_samples : int
    tolerance : float
        Fractional tolerance on bond lengths.
    min_separation : float, optional
        Minimum allowed distance between any placed site and nearby sites.
        If None, no overlap checking is performed.

    Returns
    -------
    placed_positions : np.ndarray, shape (n_sites, 3)
        Final absolute positions for all replacement sites.
    assignment : list of int
        assignment[i] = neighbor_index for active_conn_sites[i].
    """
    degree = len(active_conn_sites)

    if degree == 0:
        center = (
            np.mean(neighbor_coords, axis=0)
            if len(neighbor_coords) > 0
            else np.zeros(3)
        )
        R = _random_rotation_matrix(rng)
        placed = (R @ replacement.coordinates.T).T + center
        return placed.astype(np.float32), []

    # --- Compute target positions for connection sites ---
    neighbor_centroid = np.mean(neighbor_coords, axis=0)

    target_positions = np.zeros((degree, 3), dtype=np.float32)
    for i in range(degree):
        direction = neighbor_centroid - neighbor_coords[i]
        d = np.linalg.norm(direction)
        if d > 1e-10:
            direction_hat = direction / d
        else:
            direction_hat = rng.standard_normal(3).astype(np.float32)
            direction_hat /= np.linalg.norm(direction_hat)
        target_positions[i] = (
            neighbor_coords[i] + direction_hat * target_bond_lengths[i]
        )

    # --- Get unique connection site info ---
    unique_active = list(dict.fromkeys(active_conn_sites))
    unique_conn_positions = replacement.coordinates[unique_active]
    n_unique = len(unique_active)

    conn_to_target_indices = {}
    for i, site_idx in enumerate(active_conn_sites):
        conn_to_target_indices.setdefault(site_idx, []).append(i)

    unique_targets = np.array(
        [
            np.mean(target_positions[conn_to_target_indices[s]], axis=0)
            for s in unique_active
        ],
        dtype=np.float32,
    )

    # --- Initial Procrustes alignment ---
    if n_unique == 1:
        base_placed, base_perm, rotation_type, rotation_info = _initial_align_single(
            replacement,
            unique_active,
            unique_conn_positions,
            unique_targets,
            neighbor_centroid,
            rng,
        )
    elif n_unique >= 2:
        base_placed, base_perm, rotation_type, rotation_info = _initial_align_multi(
            replacement,
            unique_active,
            unique_conn_positions,
            unique_targets,
            rng,
        )
    else:
        R = _random_rotation_matrix(rng)
        center = np.mean(neighbor_coords, axis=0)
        placed = (R @ replacement.coordinates.T).T + center
        return placed.astype(np.float32), []

    # --- Optimize residual DOF for overlap avoidance + bond validity ---
    placed = _optimize_residual_rotation(
        base_placed,
        rotation_type,
        rotation_info,
        replacement,
        active_conn_sites,
        neighbor_coords,
        target_bond_lengths,
        nearby_coords,
        tolerance,
        min_separation,
        n_rotation_samples,
        rng,
    )

    # --- Build assignment ---
    assignment = _build_assignment_from_perm(
        active_conn_sites, unique_active, conn_to_target_indices, base_perm, degree
    )

    return placed.astype(np.float32), assignment


def _initial_align_single(
    replacement,
    unique_active,
    unique_conn_positions,
    unique_targets,
    neighbor_centroid,
    rng,
):
    """Initial alignment for 1 unique connection site.

    Returns base placement and rotation info for residual DOF.

    Returns
    -------
    base_placed : np.ndarray, shape (n_sites, 3)
    perm : list
    rotation_type : str
        'axis' for rotation about an axis, 'full' for full sphere
    rotation_info : dict
        Contains 'pivot' and 'axis' for residual rotation.
    """
    c = unique_conn_positions[0]
    c_norm = np.linalg.norm(c)

    if c_norm < 1e-10:
        R = _random_rotation_matrix(rng)
        t = unique_targets[0]
        base_placed = (R @ replacement.coordinates.T).T + t
        return base_placed, [0], "full", {"pivot": unique_targets[0]}

    c_hat = c / c_norm
    desired_dir = unique_targets[0] - neighbor_centroid
    desired_norm = np.linalg.norm(desired_dir)
    if desired_norm < 1e-10:
        desired_dir = rng.standard_normal(3).astype(np.float32)
        desired_norm = np.linalg.norm(desired_dir)
    desired_hat = desired_dir / desired_norm

    R = _rotation_between_vectors(c_hat, desired_hat)
    t = unique_targets[0] - R @ c
    base_placed = (R @ replacement.coordinates.T).T + t

    # Residual DOF: rotation about axis from center to connection site
    conn_placed = base_placed[unique_active[0]]
    center = t  # crosslinker center position
    axis = conn_placed - center
    axis_norm = np.linalg.norm(axis)
    if axis_norm > 1e-10:
        axis_hat = axis / axis_norm
    else:
        axis_hat = desired_hat

    return base_placed, [0], "axis", {"pivot": conn_placed, "axis": axis_hat}


def _initial_align_multi(
    replacement,
    unique_active,
    unique_conn_positions,
    unique_targets,
    rng,
):
    """Initial alignment for 2+ unique connection sites via Procrustes.

    Returns
    -------
    base_placed : np.ndarray
    perm : list
    rotation_type : str
        'axis' for 2 sites, 'none' for 3+
    rotation_info : dict
    """
    n_unique = len(unique_active)

    best_R = np.eye(3, dtype=np.float32)
    best_t = np.zeros(3, dtype=np.float32)
    best_cost = np.inf
    best_perm = list(range(n_unique))

    for perm in permutations(range(n_unique)):
        perm_targets = unique_targets[list(perm)]

        src_centroid = np.mean(unique_conn_positions, axis=0)
        tgt_centroid = np.mean(perm_targets, axis=0)

        src_centered = unique_conn_positions - src_centroid
        tgt_centered = perm_targets - tgt_centroid

        R = _kabsch_rotation(src_centered, tgt_centered)
        t = tgt_centroid - R @ src_centroid

        placed_conn = (R @ unique_conn_positions.T).T + t
        cost = np.sum((placed_conn - perm_targets) ** 2)

        if cost < best_cost:
            best_cost = cost
            best_R = R
            best_t = t
            best_perm = list(perm)

    base_placed = (best_R @ replacement.coordinates.T).T + best_t

    if n_unique == 2:
        # Residual DOF: rotation about axis connecting the two connection sites
        conn0 = base_placed[unique_active[0]]
        conn1 = base_placed[unique_active[1]]
        axis = conn1 - conn0
        axis_norm = np.linalg.norm(axis)
        if axis_norm > 1e-10:
            axis_hat = axis / axis_norm
        else:
            axis_hat = np.array([0, 0, 1], dtype=np.float32)
        midpoint = (conn0 + conn1) / 2.0
        return base_placed, best_perm, "axis", {"pivot": midpoint, "axis": axis_hat}
    else:
        # n_unique >= 3: fully determined, no residual DOF
        return base_placed, best_perm, "none", {}


def _optimize_residual_rotation(
    base_placed,
    rotation_type,
    rotation_info,
    replacement,
    active_conn_sites,
    neighbor_coords,
    target_bond_lengths,
    nearby_coords,
    tolerance,
    min_separation,
    n_rotation_samples,
    rng,
):
    """Optimize residual rotational DOF to avoid overlaps while maintaining bond lengths.

    Tries multiple rotations (and offsets if needed), scoring each by:
    1. Bond length validity (hard constraint)
    2. Minimum distance to nearby sites (soft maximize)

    If no rotation passes min_separation, tries pushing the crosslinker along
    the perpendicular/normal direction at fractional offsets, then re-sweeps.

    Parameters
    ----------
    base_placed : np.ndarray, shape (n_sites, 3)
        Initial placement from Procrustes.
    rotation_type : str
        'axis', 'full', or 'none'
    rotation_info : dict
        Contains 'pivot' and 'axis' as needed.
    replacement : CrosslinkerGeometry
    active_conn_sites : list of int
    neighbor_coords : np.ndarray, shape (degree, 3)
    target_bond_lengths : np.ndarray, shape (degree,)
    nearby_coords : np.ndarray, shape (M, 3)
    tolerance : float
    min_separation : float or None
    n_rotation_samples : int
    rng : numpy.random.Generator

    Returns
    -------
    best_placed : np.ndarray, shape (n_sites, 3)
    """
    if rotation_type == "none" or n_rotation_samples <= 1:
        # No residual DOF or no sampling requested
        return base_placed

    has_nearby = nearby_coords is not None and len(nearby_coords) > 0
    if not has_nearby and min_separation is None:
        return base_placed

    # Generate candidate rotations
    if rotation_type == "axis":
        pivot = rotation_info["pivot"]
        axis_hat = rotation_info["axis"]
        angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)

        def generate_candidate(angle):
            R_rot = _rotation_about_axis(axis_hat, angle)
            return (R_rot @ (base_placed - pivot).T).T + pivot

    elif rotation_type == "full":
        pivot = rotation_info["pivot"]

        # Sample random rotations about the pivot
        def generate_candidate(idx):
            R_rot = _random_rotation_matrix(rng)
            return (R_rot @ (base_placed - pivot).T).T + pivot

        angles = list(range(n_rotation_samples))
    else:
        return base_placed

    # --- Score each candidate ---
    best_placed = base_placed.copy()
    best_min_dist = -np.inf
    best_valid = False

    bond_tol = tolerance  # fractional

    for angle in angles:
        candidate = generate_candidate(angle)

        # Check bond validity
        bonds_ok = _check_bonds_valid(
            candidate, active_conn_sites, neighbor_coords, target_bond_lengths, bond_tol
        )
        if not bonds_ok:
            continue

        # Score: minimum distance to nearby sites
        if has_nearby:
            min_dist = _min_distance_to_neighbors(candidate, nearby_coords)
        else:
            min_dist = np.inf

        # Track best valid placement
        is_valid = min_separation is None or min_dist >= min_separation

        if is_valid and (not best_valid or min_dist > best_min_dist):
            best_placed = candidate.copy()
            best_min_dist = min_dist
            best_valid = True
        elif not best_valid and min_dist > best_min_dist:
            # No valid found yet, track best-effort
            best_placed = candidate.copy()
            best_min_dist = min_dist

    # --- If no valid rotation found, try offset along perpendicular ---
    if not best_valid and min_separation is not None and rotation_type == "axis":
        best_placed, best_valid = _try_offset_placements(
            base_placed,
            rotation_info,
            replacement,
            active_conn_sites,
            neighbor_coords,
            target_bond_lengths,
            nearby_coords,
            bond_tol,
            min_separation,
            n_rotation_samples,
            rng,
        )

    return best_placed


def _try_offset_placements(
    base_placed,
    rotation_info,
    replacement,
    active_conn_sites,
    neighbor_coords,
    target_bond_lengths,
    nearby_coords,
    bond_tol,
    min_separation,
    n_rotation_samples,
    rng,
):
    """Try offsetting the crosslinker perpendicular to the rotation axis.

    When the crosslinker is in the plane of the backbone and all rotations
    produce overlaps, pushing it "above" or "below" the plane can resolve
    the clash while maintaining bond lengths (connection sites rotate to
    maintain the correct distance).

    Tries offsets in both directions at increasing magnitudes.

    Parameters
    ----------
    ... (same as parent)

    Returns
    -------
    best_placed : np.ndarray, shape (n_sites, 3)
    found_valid : bool
    """
    axis_hat = rotation_info["axis"]
    pivot = rotation_info["pivot"]

    # Find a perpendicular direction to offset along
    # Use the cross product of axis with a non-parallel vector
    perp = np.array([1, 0, 0], dtype=np.float32)
    if abs(np.dot(axis_hat, perp)) > 0.9:
        perp = np.array([0, 1, 0], dtype=np.float32)
    offset_dir = np.cross(axis_hat, perp)
    offset_dir /= np.linalg.norm(offset_dir)

    # Try offsets at increasing distances
    # Max offset: limited by how much bond length can stretch within tolerance
    max_offset = np.mean(target_bond_lengths) * bond_tol * 3
    offsets = np.linspace(0.02 * max_offset, max_offset, 5)

    best_placed = base_placed.copy()
    best_min_dist = -np.inf
    found_valid = False

    for sign in [1.0, -1.0]:
        for offset_mag in offsets:
            offset = offset_dir * offset_mag * sign
            shifted = base_placed + offset

            # Re-sweep rotation at this offset
            angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)
            shifted_pivot = pivot + offset

            for angle in angles:
                R_rot = _rotation_about_axis(axis_hat, angle)
                candidate = (R_rot @ (shifted - shifted_pivot).T).T + shifted_pivot

                bonds_ok = _check_bonds_valid(
                    candidate,
                    active_conn_sites,
                    neighbor_coords,
                    target_bond_lengths,
                    bond_tol,
                )
                if not bonds_ok:
                    continue

                if nearby_coords is not None and len(nearby_coords) > 0:
                    min_dist = _min_distance_to_neighbors(candidate, nearby_coords)
                else:
                    min_dist = np.inf

                is_valid = min_dist >= min_separation

                if is_valid and (not found_valid or min_dist > best_min_dist):
                    best_placed = candidate.copy()
                    best_min_dist = min_dist
                    found_valid = True
                elif not found_valid and min_dist > best_min_dist:
                    best_placed = candidate.copy()
                    best_min_dist = min_dist

            if found_valid:
                break
        if found_valid:
            break

    return best_placed, found_valid


def _check_bonds_valid(
    placed, active_conn_sites, neighbor_coords, target_bond_lengths, tolerance
):
    """Check if all connection-to-neighbor bonds are within tolerance.

    Parameters
    ----------
    placed : np.ndarray, shape (n_sites, 3)
    active_conn_sites : list of int
    neighbor_coords : np.ndarray, shape (degree, 3)
    target_bond_lengths : np.ndarray, shape (degree,)
    tolerance : float (fractional)

    Returns
    -------
    valid : bool
    """
    for i, site_idx in enumerate(active_conn_sites):
        delta = placed[site_idx] - neighbor_coords[i]
        actual = np.linalg.norm(delta)
        target = target_bond_lengths[i]
        if abs(actual - target) > target * tolerance:
            return False
    return True


def _min_distance_to_neighbors(positions, nearby_coords):
    """Minimum distance from any position to any nearby coord."""
    if len(nearby_coords) == 0:
        return np.inf
    diffs = positions[:, None, :] - nearby_coords[None, :, :]
    dists = np.linalg.norm(diffs, axis=-1)
    return float(np.min(dists))


def _build_assignment_from_perm(
    active_conn_sites, unique_active, conn_to_target_indices, best_perm, degree
):
    """Build full assignment list from permutation over unique sites.

    Parameters
    ----------
    active_conn_sites : list of int
        Connection sites being used (length = degree).
    unique_active : list of int
        Unique site indices.
    conn_to_target_indices : dict
        site_idx -> [target indices using that site].
    best_perm : list of int
        Permutation mapping unique site position -> unique target position.
    degree : int
        Number of connections being made.

    Returns
    -------
    assignment : list of int
        assignment[i] = neighbor_index for active_conn_sites[i].
    """
    if degree == len(unique_active) and all(
        len(v) == 1 for v in conn_to_target_indices.values()
    ):
        # Simple: 1:1, no repeated sites
        # best_perm[i] maps unique site i to unique target i
        # But we need to map conn_idx -> neighbor_idx
        # active_conn_sites[conn_idx] -> unique_active.index(that) -> perm -> target
        assignment = []
        for conn_idx, site_idx in enumerate(active_conn_sites):
            unique_pos = unique_active.index(site_idx)
            assignment.append(best_perm[unique_pos])
        return assignment

    # General case with possible repeated sites
    assignment = [0] * degree
    site_counter = {}

    for conn_idx, site_idx in enumerate(active_conn_sites):
        unique_pos = unique_active.index(site_idx)
        target_group_pos = best_perm[unique_pos]
        target_site = unique_active[target_group_pos]
        target_indices = conn_to_target_indices[target_site]

        count = site_counter.get(site_idx, 0)
        assignment[conn_idx] = target_indices[count % len(target_indices)]
        site_counter[site_idx] = count + 1

    return assignment


def _verify_and_fix_bonds(
    placed_positions,
    active_conn_sites,
    assignment,
    neighbor_coords,
    target_bond_lengths,
    tolerance,
):
    """Verify connection-to-neighbor distances and nudge if outside tolerance.

    Moves each connection site along the line to its neighbor to achieve
    exactly target_bond_length if it deviates too much. Then shifts the
    entire rigid body to maintain internal geometry as much as possible.

    Parameters
    ----------
    placed_positions : np.ndarray, shape (n_sites, 3)
    active_conn_sites : list of int
    assignment : list of int
    neighbor_coords : np.ndarray, shape (degree, 3)
    target_bond_lengths : np.ndarray, shape (degree,)
    tolerance : float

    Returns
    -------
    placed_positions : np.ndarray, shape (n_sites, 3)
        Corrected positions.
    """
    placed = placed_positions.copy()
    degree = len(active_conn_sites)

    needs_correction = False
    corrections = []

    for conn_idx, neighbor_idx in enumerate(assignment):
        site_local = active_conn_sites[conn_idx]
        conn_pos = placed[site_local]
        nb_pos = neighbor_coords[neighbor_idx]
        target_bl = target_bond_lengths[conn_idx]

        delta = conn_pos - nb_pos
        actual_bl = np.linalg.norm(delta)

        if actual_bl < 1e-10:
            # Degenerate: connection site is on top of neighbor
            needs_correction = True
            corrections.append((site_local, neighbor_idx, conn_idx))
            continue

        error = abs(actual_bl - target_bl)
        if error > target_bl * tolerance:
            needs_correction = True
            corrections.append((site_local, neighbor_idx, conn_idx))

    if not needs_correction:
        return placed

    # Strategy: compute ideal connection positions and shift entire rigid body
    # to best match them (least-squares center shift)
    ideal_conn_positions = np.zeros((degree, 3), dtype=np.float32)
    for conn_idx, neighbor_idx in enumerate(assignment):
        site_local = active_conn_sites[conn_idx]
        nb_pos = neighbor_coords[neighbor_idx]
        target_bl = target_bond_lengths[conn_idx]

        conn_pos = placed[site_local]
        delta = conn_pos - nb_pos
        d = np.linalg.norm(delta)
        if d > 1e-10:
            direction = delta / d
        else:
            # Pick direction from crosslinker center to this site
            center = np.mean(placed, axis=0)
            direction = placed[site_local] - center
            d2 = np.linalg.norm(direction)
            if d2 > 1e-10:
                direction = direction / d2
            else:
                direction = np.array([1, 0, 0], dtype=np.float32)

        ideal_conn_positions[conn_idx] = nb_pos + direction * target_bl

    # Compute offset: difference between where connection sites are and should be
    current_conn_positions = np.array(
        [placed[active_conn_sites[i]] for i in range(degree)], dtype=np.float32
    )
    offsets = ideal_conn_positions - current_conn_positions
    mean_offset = np.mean(offsets, axis=0)

    # Apply rigid shift to all sites
    placed += mean_offset

    # Final per-site correction for any remaining error
    for conn_idx, neighbor_idx in enumerate(assignment):
        site_local = active_conn_sites[conn_idx]
        nb_pos = neighbor_coords[neighbor_idx]
        target_bl = target_bond_lengths[conn_idx]

        delta = placed[site_local] - nb_pos
        actual_bl = np.linalg.norm(delta)

        if actual_bl < 1e-10:
            center = np.mean(placed, axis=0)
            direction = placed[site_local] - center
            d2 = np.linalg.norm(direction)
            direction = direction / d2 if d2 > 1e-10 else np.array([1, 0, 0])
            placed[site_local] = nb_pos + direction * target_bl
        elif abs(actual_bl - target_bl) > target_bl * tolerance:
            direction = delta / actual_bl
            placed[site_local] = nb_pos + direction * target_bl

    return placed


def _get_pbc_info(volume_constraint):
    """Extract PBC flags and box lengths from volume constraint."""
    if isinstance(volume_constraint, CuboidConstraint):
        pbc = np.array(volume_constraint.pbc, dtype=bool)
        box_lengths = volume_constraint.box_lengths.astype(np.float32)
    elif isinstance(volume_constraint, CylinderConstraint):
        pbc = np.array([False, False, volume_constraint.periodic_height], dtype=bool)
        box_lengths = np.array(
            [
                volume_constraint.radius * 2,
                volume_constraint.radius * 2,
                volume_constraint.height,
            ],
            dtype=np.float32,
        )
    else:
        pbc = np.array([False, False, False], dtype=bool)
        box_lengths = np.array([np.inf, np.inf, np.inf], dtype=np.float32)
    return pbc, box_lengths


def _apply_mic(delta, pbc, box_lengths):
    """Apply minimum image convention to a displacement vector."""
    delta = np.asarray(delta, dtype=np.float32)
    for dim in range(3):
        if pbc[dim]:
            delta[dim] -= np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
    return delta


def _wrap_positions(positions, pbc, box_lengths):
    """Wrap positions into periodic box [-L/2, L/2]."""
    positions = positions.copy()
    for dim in range(3):
        if pbc[dim]:
            positions[:, dim] -= (
                np.round(positions[:, dim] / box_lengths[dim]) * box_lengths[dim]
            )
    return positions


def _get_nearby_nonbonded(
    path, site, neighbors, site_pos, radius, pbc, box_lengths, excluded_sites
):
    """Get coordinates of nearby non-bonded sites for overlap computation.

    Returns positions relative to site_pos.
    """
    neighbors_set = set(neighbors) | {site} | excluded_sites

    nearby = []
    for node in path.bond_graph.nodes():
        if node in neighbors_set:
            continue
        delta = path.coordinates[node] - site_pos
        delta = _apply_mic(delta, pbc, box_lengths)
        dist = np.linalg.norm(delta)
        if dist <= radius:
            nearby.append(delta)

    if len(nearby) == 0:
        return np.empty((0, 3), dtype=np.float32)
    return np.array(nearby, dtype=np.float32)


def crosslink(
    path,
    crosslinker=None,
    bead_name="_R",
    backbone_name="_A",
    crosslink_bond_length=0.2,
    tolerance=0.1,
    excluded_bond_depth=2,
    n_connection_sites=2,
    volume_constraint=None,
    initial_point=None,
    seed=42,
    chunk_size=512,
    run_on_gpu=False,
    n_rotation_samples=36,
    overlap_radius=None,
    min_separation=None,
):
    """
    Create a crosslink bonded to backbone beads, modifying path in place.

    Selects backbone beads whose pairwise spacing matches the crosslinker
    geometry, places a temporary placeholder, then calls ``replace_sites``
    which handles rigid-body alignment, bond length guarantees, and overlap
    avoidance in a single pass.

    Parameters
    ----------
    path : Path
        The Path object to modify in place.
    crosslinker : CrosslinkerGeometry, optional
        If None, uses single-site from ``bead_name`` and ``n_connection_sites``.
    bead_name : str, default "_R"
    backbone_name : str, default "_A"
    crosslink_bond_length : float, default 0.2
        Desired bond length from connection sites to backbone beads (nm).
    tolerance : float, default 0.1
        Fractional tolerance on bond lengths (±10%).
    excluded_bond_depth : int, default 2
    n_connection_sites : int, default 2
    volume_constraint : optional
    initial_point : int or array-like, optional
    seed : int, default 42
    chunk_size : int, default 512
    run_on_gpu : bool, default False
    n_rotation_samples : int, default 36
    overlap_radius : float, optional
    min_separation : float, optional
        Minimum distance from placed sites to existing sites.
        Defaults to ``crosslink_bond_length * 0.5``.

    Returns
    -------
    Path
        Same path, modified in place.
    """
    if crosslinker is None:
        crosslinker = CrosslinkerGeometry.single_site(
            bead_name=bead_name, n_connections=n_connection_sites
        )

    if min_separation is None:
        min_separation = crosslink_bond_length * 0.5

    n_backbone_connections = crosslinker.n_connections
    rng = np.random.default_rng(seed + len(path.coordinates))

    # --- Compute search geometry ---
    ideal_bb_positions, ideal_pairwise_distances, geom_search_radius = (
        _compute_ideal_backbone_distances(crosslinker, crosslink_bond_length)
    )

    if n_backbone_connections > 1:
        max_pair_dist = np.max(ideal_pairwise_distances)
    else:
        max_pair_dist = 0.0

    search_radius = max_pair_dist * (1 + tolerance) + crosslink_bond_length * tolerance

    # --- Find backbone targets ---
    backbone_nodes = [
        node for node in path.bond_graph.nodes() if path.beads[node] == backbone_name
    ]
    backbone_subgraph = path.bond_graph.subgraph(backbone_nodes)

    if len(backbone_nodes) == 0:
        raise ValueError(f"No backbone beads with name '{backbone_name}' found in path")
    if len(backbone_nodes) < n_backbone_connections:
        raise ValueError(
            f"Not enough backbone beads ({len(backbone_nodes)}) for "
            f"{n_backbone_connections} connection sites"
        )

    pbc, box_lengths = _get_pbc_info(volume_constraint)

    candidate_nodes = [
        node for node in backbone_nodes if path.bond_graph.degree[node] <= 2
    ]
    candidate_coords = np.array(path.coordinates[candidate_nodes], dtype=np.float32)

    def get_reference_points(path, initial_point):
        if initial_point is not None:
            if isinstance(initial_point, (int, np.integer)):
                if initial_point not in path.bond_graph.nodes:
                    raise ValueError(f"Node {initial_point} not found in bond_graph")
                return np.array([initial_point]), np.array(
                    [path.coordinates[initial_point]]
                )
            else:
                initial_point32 = np.asarray(initial_point, dtype=np.float32)
                sq_distances = calculate_sq_distances(
                    initial_point32, candidate_coords, pbc=pbc, box_lengths=box_lengths
                )
                order = np.argsort(sq_distances)
                return (
                    np.array(candidate_nodes)[order],
                    path.coordinates[np.array(candidate_nodes)[order]],
                )
        else:
            shuffled = rng.choice(
                candidate_nodes, size=len(candidate_nodes), replace=False
            )
            return shuffled, path.coordinates[shuffled]

    ref_nodes, ref_coords = get_reference_points(path, initial_point)

    found_ref = False
    selected_nodes = []

    for ref_node, ref_coord in zip(ref_nodes, ref_coords):
        selected_nodes = [ref_node]

        sq_distances = calculate_sq_distances(
            ref_coord, candidate_coords, pbc=pbc, box_lengths=box_lengths
        )
        distances = np.sqrt(sq_distances)

        within_radius_mask = distances <= search_radius
        possible_pairs = np.where(within_radius_mask)[0]

        if len(possible_pairs) < n_backbone_connections - 1:
            continue

        closest_paired_nodes = possible_pairs[np.argsort(distances[possible_pairs])]
        excluded_nodes = set(
            nx.single_source_shortest_path_length(
                backbone_subgraph, ref_node, cutoff=excluded_bond_depth
            ).keys()
        )

        for idx in closest_paired_nodes:
            node = candidate_nodes[idx]
            if node in excluded_nodes:
                continue

            is_excluded = False
            for selected in selected_nodes[1:]:
                if selected in backbone_subgraph:
                    sel_excluded = set(
                        nx.single_source_shortest_path_length(
                            backbone_subgraph, selected, cutoff=excluded_bond_depth
                        ).keys()
                    )
                    if node in sel_excluded:
                        is_excluded = True
                        break
            if is_excluded:
                continue

            selected_nodes.append(int(node))
            if len(selected_nodes) >= n_backbone_connections:
                # Validate pairwise geometry
                sel_coords = np.array(
                    [path.coordinates[n] for n in selected_nodes], dtype=np.float32
                )
                unwrapped_check = _unwrap_coords(sel_coords, pbc, box_lengths)
                valid, _ = _validate_backbone_group(
                    unwrapped_check,
                    ideal_pairwise_distances,
                    crosslink_bond_length,
                    tolerance,
                    pbc,
                    box_lengths,
                )
                if valid:
                    found_ref = True
                    break
                else:
                    selected_nodes.pop()

        if found_ref:
            break

    if not found_ref:
        n_clinks = sum(1 for b in path.beads if b in set(crosslinker.beads))
        raise PathConvergenceError(
            f"Could not find {n_backbone_connections} non-neighboring backbone beads "
            f"matching crosslinker geometry with bond_length={crosslink_bond_length} "
            f"(tolerance=±{tolerance * 100:.0f}%)."
            f"\nCurrent crosslinks: {n_clinks}. "
            "Ways to increase crosslinking:\n"
            "  - Increase crosslink_bond_length\n"
            "  - Increase tolerance\n"
            "  - Decrease excluded_bond_depth\n"
            "  - Pack at higher density\n"
            "  - Relax structure"
        )

    # --- Place placeholder at centroid ---
    selected_coords = np.array(path.coordinates[selected_nodes], dtype=np.float32)
    unwrapped = _unwrap_coords(selected_coords, pbc, box_lengths)
    crosslink_center = np.mean(unwrapped, axis=0)

    for dim in range(3):
        if pbc[dim]:
            crosslink_center[dim] -= (
                np.round(crosslink_center[dim] / box_lengths[dim]) * box_lengths[dim]
            )

    _placeholder_name = "__placeholder__"
    path.append_coordinates(crosslink_center, _placeholder_name)
    placeholder_idx = len(path.coordinates) - 1

    for backbone_node in selected_nodes:
        path.bond_graph.add_edge(
            int(placeholder_idx),
            int(backbone_node),
            bond_type=(_placeholder_name, backbone_name),
        )

    # --- Single-site shortcut ---
    if crosslinker.n_sites == 1:
        path.beads[placeholder_idx] = crosslinker.beads[0]
        for backbone_node in selected_nodes:
            path.bond_graph.edges[placeholder_idx, backbone_node]["bond_type"] = (
                str(crosslinker.beads[0]),
                backbone_name,
            )
        if "name" in path.bond_graph.nodes[placeholder_idx]:
            path.bond_graph.nodes[placeholder_idx]["name"] = str(crosslinker.beads[0])

        # Fix single-site bond length if needed
        for backbone_node in selected_nodes:
            delta = path.coordinates[placeholder_idx] - path.coordinates[backbone_node]
            delta = _apply_mic(delta, pbc, box_lengths)
            actual_bl = np.linalg.norm(delta)
            if (
                abs(actual_bl - crosslink_bond_length)
                > crosslink_bond_length * tolerance
            ):
                if actual_bl > 1e-10:
                    direction = delta / actual_bl
                else:
                    direction = rng.standard_normal(3).astype(np.float32)
                    direction /= np.linalg.norm(direction)
                path.coordinates[placeholder_idx] = path.coordinates[
                    backbone_node
                ] + _apply_mic(direction * crosslink_bond_length, pbc, box_lengths)
        return path

    # --- Multi-site: replace_sites handles everything ---
    replace_sites(
        path,
        replacement=crosslinker,
        sites=[placeholder_idx],
        bond_length=crosslink_bond_length,
        tolerance=tolerance,
        min_separation=min_separation,
        volume_constraint=volume_constraint,
        n_rotation_samples=n_rotation_samples,
        overlap_radius=overlap_radius,
        seed=seed,
    )

    return path


def _compute_ideal_backbone_distances(crosslinker, bond_length):
    """Compute ideal pairwise distances between backbone targets given geometry and bond length.

    For a crosslinker with connection sites at positions c_i (relative to center),
    the ideal backbone position for site i is:
        b_i = c_i + bond_length * (c_i / |c_i|)
    i.e., the backbone sits along the same radial direction, at distance bond_length
    beyond the connection site.

    Parameters
    ----------
    crosslinker : CrosslinkerGeometry
        The crosslinker geometry.
    bond_length : float
        Desired bond length from connection site to backbone bead.

    Returns
    -------
    ideal_bb_positions : np.ndarray, shape (n_connections, 3)
        Ideal backbone positions relative to crosslinker center.
    pairwise_distances : np.ndarray, shape (n_connections, n_connections)
        Matrix of ideal pairwise distances between backbone targets.
    search_radius : float
        Maximum distance from center to any ideal backbone position.
    """
    conn_sites = crosslinker.connection_sites
    conn_positions = crosslinker.coordinates[conn_sites]

    # Compute ideal backbone positions
    ideal_bb_positions = np.zeros_like(conn_positions)
    for i in range(len(conn_sites)):
        c_i = conn_positions[i]
        c_norm = np.linalg.norm(c_i)
        if c_norm > 1e-10:
            # Backbone sits beyond connection site along radial direction
            c_hat = c_i / c_norm
            ideal_bb_positions[i] = c_i + bond_length * c_hat
        else:
            # Connection site at origin (single-site case): backbone at bond_length in any direction
            ideal_bb_positions[i] = np.array([bond_length, 0, 0], dtype=np.float32)

    # Pairwise distances
    n = len(conn_sites)
    pairwise_distances = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(ideal_bb_positions[i] - ideal_bb_positions[j])
            pairwise_distances[i, j] = d
            pairwise_distances[j, i] = d

    # Search radius: max distance from center to any backbone
    search_radius = np.max(np.linalg.norm(ideal_bb_positions, axis=1))

    return ideal_bb_positions, pairwise_distances, search_radius


def _unwrap_coords(coords, pbc, box_lengths):
    """Unwrap coordinates relative to the first point."""
    unwrapped = np.empty_like(coords)
    unwrapped[0] = coords[0]
    ref = coords[0]
    for k in range(1, len(coords)):
        delta = coords[k] - ref
        for dim in range(3):
            if pbc[dim]:
                delta[dim] -= np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
        unwrapped[k] = ref + delta
    return unwrapped


def _compute_optimal_center(
    crosslinker, backbone_coords, bond_length, pbc, box_lengths
):
    """Compute the optimal crosslinker center so connection sites are at bond_length from backbones.

    The center is placed such that each connection site (after rotation) sits
    at exactly bond_length from its target backbone bead.

    For the ideal placement:
        backbone_i = center + c_i_hat * (|c_i| + bond_length)

    So: center = backbone_i - c_i_hat * (|c_i| + bond_length)

    We average over all connection sites for robustness.

    Parameters
    ----------
    crosslinker : CrosslinkerGeometry
        The crosslinker.
    backbone_coords : np.ndarray, shape (n_connections, 3)
        Unwrapped backbone coordinates.
    bond_length : float
        Desired bond length.
    pbc, box_lengths : PBC info.

    Returns
    -------
    center : np.ndarray, shape (3,)
    """
    conn_sites = crosslinker.connection_sites
    conn_positions = crosslinker.coordinates[conn_sites]
    n = len(conn_sites)

    if n == 0:
        return np.mean(backbone_coords, axis=0)

    # For single-site crosslinker (all connection sites at same position):
    if np.allclose(conn_positions, conn_positions[0], atol=1e-8):
        # Center is at bond_length distance from centroid toward... just use centroid
        return np.mean(backbone_coords, axis=0)

    # Compute directions from crosslinker center to each connection site
    conn_norms = np.linalg.norm(conn_positions, axis=1, keepdims=True)
    conn_hat = conn_positions / np.maximum(conn_norms, 1e-10)

    # Ideal: backbone_i is at direction c_i_hat, distance |c_i| + bond_length from center
    # So center = backbone_i - c_i_hat * (|c_i| + bond_length)
    # But we haven't rotated the crosslinker yet, so we need to figure out the rotation.

    # Target vectors: from centroid of backbones to each backbone
    bb_centroid = np.mean(backbone_coords, axis=0)
    target_vectors = backbone_coords - bb_centroid
    target_norms = np.linalg.norm(target_vectors, axis=1, keepdims=True)
    target_hat = target_vectors / np.maximum(target_norms, 1e-10)

    # The ideal distance from center to each backbone:
    # |c_i| + bond_length
    ideal_radii = conn_norms.flatten() + bond_length

    # Center = backbone_i - ideal_radius_i * target_hat_i (if we orient along target directions)
    # Average:
    centers = backbone_coords - ideal_radii[:, None] * target_hat
    center = np.mean(centers, axis=0)

    return center.astype(np.float32)


def _validate_backbone_group(
    selected_coords,
    ideal_pairwise_distances,
    bond_length,
    tolerance,
    pbc,
    box_lengths,
):
    """Check if backbone beads match expected geometry within tolerance.

    Returns
    -------
    valid : bool
    best_perm : list of int or None
    """
    n = len(selected_coords)
    if n <= 1:
        return True, [0] if n == 1 else []

    abs_tol = bond_length * tolerance

    # Compute actual pairwise distances (PBC-aware)
    actual_distances = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            delta = selected_coords[j] - selected_coords[i]
            for dim in range(3):
                if pbc[dim]:
                    delta[dim] -= (
                        np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
                    )
            d = np.linalg.norm(delta)
            actual_distances[i, j] = d
            actual_distances[j, i] = d

    best_perm = None
    best_cost = np.inf

    for perm in permutations(range(n)):
        valid = True
        cost = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                ideal_d = ideal_pairwise_distances[i, j]
                actual_d = actual_distances[perm[i], perm[j]]
                error = abs(actual_d - ideal_d)
                if error > abs_tol:
                    valid = False
                    break
                cost += error
            if not valid:
                break

        if valid and cost < best_cost:
            best_cost = cost
            best_perm = list(perm)

    return best_perm is not None, best_perm
