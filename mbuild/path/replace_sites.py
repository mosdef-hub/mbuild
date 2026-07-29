from itertools import permutations

import networkx as nx
import numpy as np

from mbuild.path.constraints import CuboidConstraint, CylinderConstraint


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


def _apply_mic(delta, pbc, box_lengths):
    """Apply minimum image convention to a displacement vector."""
    delta = np.asarray(delta, dtype=np.float32)
    for dim in range(3):
        if pbc[dim]:
            delta[dim] -= np.round(delta[dim] / box_lengths[dim]) * box_lengths[dim]
    return delta


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


def _wrap_positions(positions, pbc, box_lengths):
    """Wrap positions into periodic box [-L/2, L/2]."""
    positions = positions.copy()
    for dim in range(3):
        if pbc[dim]:
            positions[:, dim] -= (
                np.round(positions[:, dim] / box_lengths[dim]) * box_lengths[dim]
            )
    return positions


def _random_rotation_matrix(rng):
    """Uniform random rotation via QR decomposition."""
    H = rng.standard_normal((3, 3)).astype(np.float32)
    Q, R = np.linalg.qr(H)
    Q *= np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q.astype(np.float32)


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


def _kabsch_rotation(P, Q):
    """Optimal rotation aligning P onto Q (Kabsch algorithm)."""
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, d])
    R = Vt.T @ sign_matrix @ U.T
    return R.astype(np.float32)


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


def _get_perpendicular(axis):
    """Get a unit vector perpendicular to axis."""
    axis = np.asarray(axis, dtype=np.float32)
    perp = np.array([1, 0, 0], dtype=np.float32)
    if abs(np.dot(axis, perp)) > 0.9:
        perp = np.array([0, 1, 0], dtype=np.float32)
    perp = perp - np.dot(perp, axis) * axis
    perp /= np.linalg.norm(perp)
    return perp
