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
    ...
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
    if isinstance(sites, int):
        sites = [sites]
    else:
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
    sites_set = set(sites)

    # --- Compute all replacements ---
    replacement_data = {}

    for site in sites:
        neighbors = list(path.bond_graph.neighbors(site))
        site_pos = path.coordinates[site]
        degree = len(neighbors)

        # Unwrapped neighbor positions
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

        active_conn_sites = replacement.connection_sites[:degree]

        # Indices that are bonded to this site (excluded from overlap check)
        bonded_indices = set(neighbors) | {site} | sites_set

        # Place via rigid body alignment + overlap optimization
        placed_positions, assignment = _place_replacement_rigid_full(
            replacement,
            active_conn_sites,
            neighbor_coords,
            target_bond_lengths,
            path.coordinates,  # pass ALL coordinates
            bonded_indices,
            pbc,
            box_lengths,
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

    # --- Rebuild path in place (same as before) ---
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

    new_edges = []

    for u, v, data in path.bond_graph.edges(data=True):
        if u in sites_set or v in sites_set:
            continue
        new_edges.append((old_to_new[u], old_to_new[v], data))

    for site in sites:
        indices = site_new_indices[site]
        for i, j in replacement.bond_graph.edges():
            edge_data = {
                "bond_type": (str(replacement.beads[i]), str(replacement.beads[j]))
            }
            new_edges.append((indices[i], indices[j], edge_data))

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
                )
            }
            new_edges.append((repl_new_node, ext_new_node, edge_data))

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


def crosslink(
    path,
    crosslinker=None,
    bead_name="_R",
    backbone_name="_A",
    crosslink_bond_length=0.2,
    max_backbone_degree=2,
    tolerance=0.1,
    excluded_bond_depth=2,
    n_connection_sites=2,
    volume_constraint=None,
    initial_point=None,
    seed=42,
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
        min_separation = crosslink_bond_length * 0.2

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

    search_radius = max_pair_dist * (1 + tolerance)

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
        node
        for node in backbone_nodes
        if path.bond_graph.degree[node] <= max_backbone_degree
    ]
    candidate_coords = np.array(path.coordinates[candidate_nodes], dtype=np.float32)

    ref_nodes, ref_coords = _get_reference_points(
        path, initial_point, candidate_nodes, pbc, box_lengths, rng
    )

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
                valid_overlaps, clink_centroid = _validate_clinker(
                    path,
                    crosslinker,
                    crosslink_bond_length,
                    selected_nodes,
                    overlap_radius=min_separation,
                )
                if valid and valid_overlaps:
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

    for dim in range(3):
        if pbc[dim]:
            clink_centroid[dim] -= (
                np.round(clink_centroid[dim] / box_lengths[dim]) * box_lengths[dim]
            )

    _placeholder_name = "__placeholder__"
    path.append_coordinates(clink_centroid, _placeholder_name)
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


def _get_reference_points(path, initial_point, candidate_nodes, pbc, box_lengths, rng):
    candidate_coords = np.array(path.coordinates[candidate_nodes], dtype=np.float32)
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
        shuffled = rng.choice(candidate_nodes, size=len(candidate_nodes), replace=False)
        return shuffled, path.coordinates[shuffled]


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

    # Maximum valid pairwise distances (geometric feasibility bound)
    n = len(conn_sites)
    pairwise_distances = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            d_internal = np.linalg.norm(conn_positions[i] - conn_positions[j])
            max_d = d_internal + 2 * bond_length
            pairwise_distances[i, j] = max_d
            pairwise_distances[j, i] = max_d

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
                max_d = ideal_pairwise_distances[i, j] * (1 + tolerance)
                actual_d = actual_distances[perm[i], perm[j]]
                if actual_d > max_d:
                    valid = False
                    break
                cost += actual_d
            if not valid:
                break

        if valid and cost < best_cost:
            best_cost = cost
            best_perm = list(perm)

    return best_perm is not None, best_perm


def _validate_clinker(
    path,
    crosslinker,
    crosslink_bond_length,
    selected_nodes,
    overlap_radius=0.1,
    n_samples=200,
    pbc=None,
    box_lengths=None,
    seed=None,
    tolerance=1e-5,
):
    """Verify a point exists within constraints to place the centroid of the crosslinker.

    Treats the crosslinker as a cylindrical bead whose radius is the maximum
    distance from its centroid to any constituent site. Samples candidate
    center positions within the geometric feasibility region (intersection of
    spheres around selected backbone nodes) and uses ``check_path`` to verify
    no hard-sphere overlaps exist with the rest of the path.

    Parameters
    ----------
    path : Path
        The current path with all existing coordinates.
    crosslinker : CrosslinkerGeometry
        The crosslinker geometry to be placed.
    crosslink_bond_length : float
        Desired bond length from crosslinker connection sites to backbone beads.
    selected_nodes : list of int
        Backbone node indices that the crosslinker would connect to.
    overlap_radius : float, default 0.1
        Additional exclusion radius beyond the crosslinker's own bead radius.
        The total check radius passed to ``check_path`` is
        ``bead_radius + overlap_radius``.
    n_samples : int, default 200
        Number of candidate center positions to test.
    pbc : np.ndarray of bool or None
        Periodic boundary conditions per axis. If None, assumes no PBC.
    box_lengths : np.ndarray or None
        Box dimensions for minimum image convention.
    seed : int or None
        Random seed for reproducibility.
    tolerance : float, default 1e-5
        Tolerance passed to ``check_path`` for rounding in distance checks.

    Returns
    -------
    valid : bool
        True if at least one overlap-free center position exists.
    best_center : np.ndarray, shape (3,) or None
        The first valid candidate center position found, or the best
        (maximum clearance) candidate if none fully pass.
    """
    from mbuild.path.path_utils import check_path

    if pbc is None:
        pbc = np.array([False, False, False])
    if box_lengths is None:
        box_lengths = np.array([0.0, 0.0, 0.0], dtype=np.float32)

    rng = np.random.default_rng(seed)

    # --- Compute crosslinker bead radius (cylindrical approximation) ---
    crosslinker_centroid = np.mean(crosslinker.coordinates, axis=0)
    site_offsets = crosslinker.coordinates - crosslinker_centroid
    bead_radius = float(np.max(np.linalg.norm(site_offsets, axis=1)))

    # --- Max distance from centroid to any connection site ---
    conn_sites = crosslinker.connection_sites
    conn_positions = crosslinker.coordinates[conn_sites]
    conn_offsets = conn_positions - crosslinker_centroid
    max_conn_radius = float(np.max(np.linalg.norm(conn_offsets, axis=1)))

    # --- Feasibility: centroid must be within this distance of each backbone node ---
    feasibility_radius = crosslink_bond_length + max_conn_radius

    # --- Total hard-sphere check radius ---
    check_radius = bead_radius + overlap_radius

    # --- Build existing_points array (all path coords except selected nodes) ---
    n_total = len(path.coordinates)
    selected_set = set(selected_nodes)
    mask = np.ones(n_total, dtype=bool)
    for node in selected_set:
        if 0 <= node < n_total:
            mask[node] = False
    existing_points = np.ascontiguousarray(path.coordinates[mask], dtype=np.float32)
    # existing_points = path.coordinates

    # --- Unwrap selected node coordinates relative to first ---
    selected_coords = np.array(
        [path.coordinates[n] for n in selected_nodes], dtype=np.float32
    )
    if len(selected_coords) > 1:
        ref = selected_coords[0].copy()
        for k in range(1, len(selected_coords)):
            delta = selected_coords[k] - ref
            for dim in range(3):
                if pbc[dim]:
                    delta[dim] -= (
                        np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
                    )
            selected_coords[k] = ref + delta

    centroid_of_selected = np.mean(selected_coords, axis=0)

    # --- Compute sampling sphere radius ---
    max_dist_to_centroid = float(
        np.max(np.linalg.norm(selected_coords - centroid_of_selected, axis=1))
    )
    sample_radius = feasibility_radius - max_dist_to_centroid
    if sample_radius <= 0:
        sample_radius = feasibility_radius * 0.05

    # --- Generate candidate center positions ---
    candidates = np.zeros((n_samples, 3), dtype=np.float32)
    candidates[0] = centroid_of_selected
    for i in range(1, n_samples):
        r = sample_radius * (rng.random() ** (1.0 / 3.0))
        direction = rng.standard_normal(3).astype(np.float32)
        norm = np.linalg.norm(direction)
        if norm > 1e-10:
            direction /= norm
        else:
            direction = np.array([1.0, 0.0, 0.0], dtype=np.float32)
        candidates[i] = centroid_of_selected + r * direction

    # --- Evaluate each candidate ---
    best_center = None

    for cand in candidates:
        # Check feasibility: centroid within reach of all selected backbone nodes
        feasible = True
        for sc in selected_coords:
            delta = cand - sc
            for dim in range(3):
                if pbc[dim]:
                    delta[dim] -= (
                        np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
                    )
            dist = np.linalg.norm(delta)
            if dist > feasibility_radius:
                feasible = False
                break
        if not feasible:
            continue

        # Use check_path for hard-sphere overlap detection
        new_point = np.ascontiguousarray(cand, dtype=np.float32)
        no_overlap = check_path(
            existing_points=existing_points,
            new_point=new_point,
            radius=check_radius,
            tolerance=tolerance,
        )

        if no_overlap:
            return True, cand

        # Track first feasible candidate as fallback
        if best_center is None:
            best_center = cand.copy()

    return False, best_center


def _place_replacement_rigid_full(
    replacement,
    active_conn_sites,
    neighbor_coords,
    target_bond_lengths,
    all_path_coords,
    bonded_indices,
    pbc,
    box_lengths,
    rng,
    n_rotation_samples,
    tolerance,
    min_separation=None,
):
    """Place replacement with full-system overlap checking.

    Unlike the radius-limited version, this checks placed sites against ALL
    existing coordinates in the path (minus bonded neighbors).

    Parameters
    ----------
    replacement : CrosslinkerGeometry
    active_conn_sites : list of int
    neighbor_coords : np.ndarray, shape (degree, 3)
    target_bond_lengths : np.ndarray, shape (degree,)
    all_path_coords : np.ndarray, shape (N, 3)
        ALL coordinates in the current path.
    bonded_indices : set of int
        Indices to exclude from overlap checking.
    pbc : array-like of bool
    box_lengths : np.ndarray
    rng : numpy.random.Generator
    n_rotation_samples : int
    tolerance : float
    min_separation : float or None

    Returns
    -------
    placed_positions : np.ndarray, shape (n_sites, 3)
    assignment : list of int
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

    # --- Compute target positions ---
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

    # --- Initial alignment ---
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
        placed = (R @ replacement.coordinates.T).T + neighbor_centroid
        assignment = _build_assignment_from_perm(
            active_conn_sites, unique_active, conn_to_target_indices, [0], degree
        )
        return placed.astype(np.float32), assignment

    # --- Precompute check coords (all minus bonded) ---
    n_existing = len(all_path_coords)
    check_mask = np.ones(n_existing, dtype=bool)
    for idx in bonded_indices:
        if 0 <= idx < n_existing:
            check_mask[idx] = False
    check_coords = all_path_coords[check_mask]

    # --- Sweep rotations with full overlap check ---
    best_placed, best_min_dist, found_valid = _sweep_rotations_full(
        base_placed,
        rotation_type,
        rotation_info,
        active_conn_sites,
        neighbor_coords,
        target_bond_lengths,
        check_coords,
        pbc,
        box_lengths,
        tolerance,
        min_separation,
        n_rotation_samples,
        rng,
    )

    # --- If not valid, try offsets ---
    if not found_valid and min_separation is not None and rotation_type == "axis":
        best_placed, found_valid = _try_offsets_full(
            base_placed,
            rotation_info,
            replacement,
            active_conn_sites,
            neighbor_coords,
            target_bond_lengths,
            check_coords,
            pbc,
            box_lengths,
            tolerance,
            min_separation,
            n_rotation_samples,
            rng,
        )

    # Build assignment
    assignment = _build_assignment_from_perm(
        active_conn_sites, unique_active, conn_to_target_indices, base_perm, degree
    )

    return best_placed.astype(np.float32), assignment


def _sweep_rotations_full(
    base_placed,
    rotation_type,
    rotation_info,
    active_conn_sites,
    neighbor_coords,
    target_bond_lengths,
    check_coords,
    pbc,
    box_lengths,
    tolerance,
    min_separation,
    n_rotation_samples,
    rng,
):
    """Sweep rotations checking against full coordinate set.

    Returns
    -------
    best_placed : np.ndarray
    best_min_dist : float
    found_valid : bool
    """
    best_placed = base_placed.copy()
    best_min_dist = _min_dist_pbc(base_placed, check_coords, pbc, box_lengths)
    found_valid = min_separation is None or best_min_dist >= min_separation

    # Check if base placement bonds are valid
    if not _check_bonds_valid(
        base_placed, active_conn_sites, neighbor_coords, target_bond_lengths, tolerance
    ):
        best_min_dist = -np.inf
        found_valid = False

    if rotation_type == "axis":
        pivot = rotation_info["pivot"]
        axis_hat = rotation_info["axis"]
        angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)

        for angle in angles:
            R_rot = _rotation_about_axis(axis_hat, angle)
            candidate = (R_rot @ (base_placed - pivot).T).T + pivot

            if not _check_bonds_valid(
                candidate,
                active_conn_sites,
                neighbor_coords,
                target_bond_lengths,
                tolerance,
            ):
                continue

            min_dist = _min_dist_pbc(candidate, check_coords, pbc, box_lengths)

            if min_separation is not None and min_dist >= min_separation:
                if not found_valid or min_dist > best_min_dist:
                    best_placed = candidate.copy()
                    best_min_dist = min_dist
                    found_valid = True
            elif not found_valid and min_dist > best_min_dist:
                best_placed = candidate.copy()
                best_min_dist = min_dist

    elif rotation_type == "full":
        pivot = rotation_info.get("pivot", np.mean(base_placed, axis=0))
        for _ in range(n_rotation_samples):
            R = _random_rotation_matrix(rng)
            candidate = (R @ (base_placed - pivot).T).T + pivot

            if not _check_bonds_valid(
                candidate,
                active_conn_sites,
                neighbor_coords,
                target_bond_lengths,
                tolerance,
            ):
                continue

            min_dist = _min_dist_pbc(candidate, check_coords, pbc, box_lengths)

            if min_separation is not None and min_dist >= min_separation:
                if not found_valid or min_dist > best_min_dist:
                    best_placed = candidate.copy()
                    best_min_dist = min_dist
                    found_valid = True
            elif not found_valid and min_dist > best_min_dist:
                best_placed = candidate.copy()
                best_min_dist = min_dist

    return best_placed, best_min_dist, found_valid


def _try_offsets_full(
    base_placed,
    rotation_info,
    replacement,
    active_conn_sites,
    neighbor_coords,
    target_bond_lengths,
    check_coords,
    pbc,
    box_lengths,
    tolerance,
    min_separation,
    n_rotation_samples,
    rng,
):
    """Try offset placements with full overlap checking."""
    axis_hat = rotation_info["axis"]
    pivot = rotation_info["pivot"]

    perp1 = _get_perpendicular(axis_hat)
    perp2 = np.cross(axis_hat, perp1).astype(np.float32)
    perp2 /= np.linalg.norm(perp2)

    max_extent = float(np.max(np.linalg.norm(replacement.coordinates, axis=1)))
    offset_magnitudes = np.linspace(0.05, max_extent * 1.5, 6)

    best_placed = base_placed.copy()
    best_min_dist = _min_dist_pbc(base_placed, check_coords, pbc, box_lengths)
    found_valid = False

    for offset_mag in offset_magnitudes:
        for sign in [1.0, -1.0]:
            for perp_dir in [perp1, perp2]:
                offset = perp_dir * offset_mag * sign

                # Shift the whole placement
                shifted = base_placed + offset
                shifted_pivot = pivot + offset

                # Sweep rotations at this offset
                angles = np.linspace(0, 2 * np.pi, n_rotation_samples, endpoint=False)
                for angle in angles:
                    R_rot = _rotation_about_axis(axis_hat, angle)
                    candidate = (R_rot @ (shifted - shifted_pivot).T).T + shifted_pivot

                    if not _check_bonds_valid(
                        candidate,
                        active_conn_sites,
                        neighbor_coords,
                        target_bond_lengths,
                        tolerance,
                    ):
                        continue

                    min_dist = _min_dist_pbc(candidate, check_coords, pbc, box_lengths)

                    if min_dist >= min_separation:
                        if not found_valid or min_dist > best_min_dist:
                            best_placed = candidate.copy()
                            best_min_dist = min_dist
                            found_valid = True
                    elif not found_valid and min_dist > best_min_dist:
                        best_placed = candidate.copy()
                        best_min_dist = min_dist

                if found_valid:
                    return best_placed, found_valid

    return best_placed, found_valid


def _min_dist_pbc(positions, check_coords, pbc, box_lengths):
    """Minimum distance from any position to any check coord, with PBC."""
    if len(check_coords) == 0:
        return np.inf

    global_min_sq = np.inf
    for site_pos in positions:
        deltas = check_coords - site_pos
        for dim in range(3):
            if pbc[dim]:
                deltas[:, dim] -= (
                    np.round(deltas[:, dim] / box_lengths[dim]) * box_lengths[dim]
                )
        sq_dists = np.sum(deltas * deltas, axis=1)
        local_min = float(np.min(sq_dists))
        if local_min < global_min_sq:
            global_min_sq = local_min

    return float(np.sqrt(global_min_sq))


def _get_perpendicular(axis):
    """Get a unit vector perpendicular to axis."""
    axis = np.asarray(axis, dtype=np.float32)
    perp = np.array([1, 0, 0], dtype=np.float32)
    if abs(np.dot(axis, perp)) > 0.9:
        perp = np.array([0, 1, 0], dtype=np.float32)
    perp = perp - np.dot(perp, axis) * axis
    perp /= np.linalg.norm(perp)
    return perp
