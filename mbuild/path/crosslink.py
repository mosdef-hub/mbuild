"""Crosslinking module for polymer path structures.

Handles placement of crosslinker geometries between backbone beads,
with proper overlap avoidance and bond length enforcement under PBC.
"""

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from mbuild.exceptions import PathConvergenceError
from mbuild.path.constraints import CuboidConstraint
from mbuild.path.path_utils import check_path


class CrosslinkerGeometry:
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
    isrigid : bool
        If the crosslinker should be held rigid
    """

    def __init__(
        self,
        coordinates,
        bond_graph=None,
        bead_name="_R",
        connection_sites=None,
        isrigid=False,
    ):
        coordinates = np.asarray(coordinates, dtype=np.float32)
        self.coordinates = coordinates - coordinates.mean(axis=0)
        self.isrigid = isrigid

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
        elif isinstance(connection_sites, int):
            self.connection_sites = [connection_sites, connection_sites]
        else:
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


def crosslink(
    path,
    crosslinker=None,
    n_crosslinks=1,
    crosslink_bond_length=0.2,
    max_backbone_degree=4,
    tolerance=0.1,
    minimum_intrachain_depth=2,
    volume_constraint=None,
    minimum_separation=0.1,
    backbone_bead_name=None,
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
    if volume_constraint is None:
        box = CuboidConstraint(np.inf, np.inf, np.inf)  # aperiodic inf box is default
    elif isinstance(volume_constraint, CuboidConstraint):
        box = volume_constraint
    else:
        raise ValueError(
            f"Only supporting crosslinking in CuboidConstrint for the volume_constraint, not {volume_constraint}"
        )

    # filter down to valid backbone sites here
    backbone_sitesArray = _filter_initial_backbones(
        path, backbone_bead_name, max_backbone_degree
    )

    # key to access supported crosslink geometry functions
    if crosslinker is None:
        crosslinkKey = 0
    else:
        # use a set to account for duplicates
        crosslinkKey = len(
            set(crosslinker.connection_sites)
        )  # i.e. [0,0] means only one unique bonding site
    if isinstance(crosslink_bond_length, (float, int)) and crosslinker is not None:
        crosslink_bond_lengthList = [crosslink_bond_length] * len(
            crosslinker.connection_sites
        )
    elif isinstance(crosslink_bond_length, (float, int)):
        crosslink_bond_lengthList = [crosslink_bond_length]
    elif isinstance(crosslink_bond_length, (list, tuple, np.ndarray)):
        crosslink_bond_lengthList = list(crosslink_bond_lengthList)
        if not len(crosslink_bond_lengthList) == len(crosslinker.connection_sites):
            raise PathConvergenceError(
                f"backbone_bond_lengths list doesn't match crosslinker.connection_sites"
                f"{crosslink_bond_length}"
                and {crosslinker.connection_sites}
            )

    rng = np.random.default_rng(seed)
    n_placed = 0
    max_attempts = n_crosslinks * 100
    attempt = 0

    while n_placed < n_crosslinks and attempt < max_attempts:
        attempt += 1

        if crosslinkKey == 0:
            candidatesGenerator = direct_crosslink
            generator_kwargs = {
                "coords": path.coordinates,
                "backbone_sites": backbone_sitesArray,
                "crosslink_bond_lengthList": crosslink_bond_lengthList,
                "bond_tolerance": tolerance,
                "rng": rng,
                "bond_graph": path.bond_graph,
                "minimum_intrachain_depth": minimum_intrachain_depth,
                "box": box,
            }
        elif crosslinkKey == 1:
            candidatesGenerator = one_site_crosslink
            generator_kwargs = {
                "coords": path.coordinates,
                "backbone_sites": backbone_sitesArray,
                "crosslink_bond_lengthList": crosslink_bond_lengthList,
                "bond_tolerance": tolerance,
                "rng": rng,
                "bond_graph": path.bond_graph,
                "minimum_intrachain_depth": minimum_intrachain_depth,
                "box": box,
            }
        elif crosslinkKey == 2:
            raise ValueError(f"Unsupported for {crosslinker} as of now.")
        else:
            candidatesGenerator = None  # 3 sites and more.

        # validate generated sites
        crosslinker_pos = None  # stays None if we never generate viable candidates
        for candidate in candidatesGenerator(**generator_kwargs):
            if crosslinkKey == 0:
                crosslinker_pos = _verify_direct_bond(
                    candidate,
                    path.coordinates,
                    crosslink_bond_lengthList[0],
                    tolerance,
                    box,
                )
            elif crosslinkKey == 1:
                crosslinker_pos = _verify_one_site_bond(
                    candidate,
                    crosslinker,
                    path.coordinates,
                    minimum_separation,
                    crosslink_bond_lengthList,
                    box=box,
                    tolerance=tolerance,
                    rng=rng,
                    n_position_samples=8,  # TODO: expose to top level api
                    n_rotation_samples=12,
                )
            if crosslinker_pos is not False:
                break
        if crosslinker_pos is False:  # found a valid site, but couldn't fit overlaps
            raise PathConvergenceError(
                f"Could not find crosslinker placement with "
                f"bond_length={crosslink_bond_length} (tolerance=±{tolerance * 100:.0f}%).\n"
                f"minimum_separation={minimum_separation}.\n"
                f"minimum_intrachain_depth={minimum_intrachain_depth}.\n"
                f"Ended after {n_placed} crosslinks placed.\n"
                "Ways to increase crosslinking:\n"
                " - Modify crosslink_bond_length or tolerance\n"
                " - Increase max_backbone_degree\n"
                " - Decrease excluded_bond_depth\n"
                " - Pack at different density"
                " - Run relaxation after this error and try again/"
            )
        elif crosslinker_pos is None:  # never found a viable candidate after loop
            raise PathConvergenceError(
                "Found no valid candidate sites in backbone for crosslinking."
            )
        if crosslinkKey == 0:
            path.bond_graph.add_edge(*candidate)
        else:  # any other number of connection sites
            _insert_crosslinker(path, crosslinker, candidate, crosslinker_pos, box)
        n_placed += 1
        if verbose and n_placed % 10 == 0:
            print(f"Placed {n_placed}/{n_crosslinks} crosslinks")

    if n_placed < n_crosslinks:
        raise PathConvergenceError(
            f"Could only place {n_placed}/{n_crosslinks} crosslinks after "
            f"{attempt} attempts.\n"
            f"Consider increasing tolerance, n_rotation_samples, or "
            f"decreasing minimum_separation."
        )

    return path


# =============================================================================
# General crosslink functions
# =============================================================================


def _filter_initial_backbones(path, backbone_name, backbone_degree):
    """Return a list of valid backbones in path.

    Parameters
    ----------
    path : mbuild.path.build.Path object
        Path to use for searching backbones
    backbone_name : str or list or set
        names filter to get indices of valid backbone sites
    backbone_degree : int or tuple or list
        Specification of site backbone degree.
        If is an integer, set to the maximum allowed degree.
        If an iterable is passed, assume the first value is the minumimum
        and the second is the maximum value, inclusive.

    Returns
    -------
    backbones : np.ndarray (M,3) of integers
        Where M is the number of valid backbones. Each integer is an index in the overall path
        that can be considered a backbone.

    Notes
    -----
    Handles a single bead name, and the maximum_degree. In the future might more specifically filter
    out other features of backbones.
    """
    if backbone_name is None:
        return np.arange(len(path))
    if isinstance(backbone_name, str):
        backbone_namesSet = set(
            (backbone_name,)
        )  # every connection bond to same crosslink
    else:
        backbone_namesSet = set(backbone_name)  # assume it's some iterable
    if isinstance(backbone_degree, int):
        min_backbone_degree = 1
        max_backbone_degree = backbone_degree
    elif isinstance(backbone_degree, (tuple, list, np.ndarray)):
        min_backbone_degree = backbone_degree[0]
        max_backbone_degree = backbone_degree[1]
    n = path.bond_graph.number_of_nodes()
    degrees = np.zeros(n, dtype=int)
    nodes, degs = zip(*path.bond_graph.degree())
    degrees[list(nodes)] = degs

    name_ok = np.isin(
        path.beads, np.fromiter(backbone_namesSet, dtype=path.beads.dtype)
    )
    return np.flatnonzero(
        (degrees >= min_backbone_degree) & (degrees <= max_backbone_degree) & name_ok
    )


def _insert_crosslinker(
    path,
    crosslinker,
    candidate,
    crosslinker_pos,
    box,
):
    """Modify path inplace with crosslinker positions, and update backbone info for next cycle.

    Returns
    -------
    None
    """
    box_center = box.center

    # Wrap final positions into the box one more time as a safety net.
    if any(box.pbc):
        mask = box.pbc & ~np.isinf(box.box_lengths)
        if np.any(mask):
            shifted = (crosslinker_pos - box_center).astype(np.float32)
            wrapped = box.wrap(shifted)  # TODO: Add a cuboid wrap function
            # Apply mask to dimensions (axis 1), not particles (axis 0)
            crosslinker_pos[:, mask] = (wrapped[:, mask] + box_center[mask]).astype(
                np.float32
            )

    previous_n_nodes = len(path)
    path.add_path(crosslinker, crosslinker_pos)
    for clink_site, backbone_site in zip(crosslinker.connection_sites, candidate):
        path.bond_graph.add_edge(previous_n_nodes + clink_site, backbone_site)


# =============================================================================
# Direct Crosslinking (backbone to backbone)
# =============================================================================


def direct_crosslink(
    coords,
    backbone_sites,
    crosslink_bond_lengthList,
    bond_tolerance,
    rng,
    box,
    bond_graph=None,
    minimum_intrachain_depth=2,
):
    """Yield candidate backbone site pairs for a direct crosslink bond."""
    box_center = box.center

    bond_length = crosslink_bond_lengthList[0]
    r_max = bond_length * (1 + bond_tolerance) + 1e-4  # 1e-4 edge padding
    r_min = bond_length * (1 - bond_tolerance) - 1e-4  # edge padding
    n_backbones = len(backbone_sites)

    # Boolean mask so we can intersect the distance window with valid sites.
    backbone_set = set(backbone_sites.tolist())

    for site in rng.choice(np.arange(n_backbones), size=n_backbones, replace=False):
        origin = backbone_sites[site]

        if bond_graph is not None:
            excluded = set(
                nx.single_source_shortest_path_length(
                    bond_graph, origin, cutoff=minimum_intrachain_depth
                ).keys()
            )
        else:
            excluded = {origin}

        # _find_neighbors dispatches to AABBQuery (all-periodic) or tiled
        # cKDTree (any aperiodic axis).
        neighbors = _find_neighbors(
            coords - box_center,
            box,
            origin,
            r_max=r_max,
            r_min=r_min,
        )

        for candidate_site in neighbors:
            if candidate_site not in backbone_set or candidate_site in excluded:
                continue
            yield (origin, int(candidate_site))


def _verify_direct_bond(
    candidate,
    coords,
    bond_length,
    tolerance,
    box,
):
    """Verify a candidate backbone pair is within tolerance of the target
    direct bond length, using minimum-image distance under PBC.
    """
    pbc = box.pbc
    box_lengths = box.box_lengths

    p1, p2 = coords[list(candidate)]
    d = np.linalg.norm(_minimum_image(p2 - p1, pbc, box_lengths))

    if bond_length * (1 - tolerance) <= d <= bond_length * (1 + tolerance):
        return d


# =============================================================================
# One-Site Crosslinker
# =============================================================================


def one_site_crosslink(
    coords,
    backbone_sites,
    crosslink_bond_lengthList,
    bond_tolerance,
    rng,
    box,
    bond_graph=None,
    minimum_intrachain_depth=2,
):
    """Yield candidate backbone site pairs for a one-site crosslink."""
    box_center = box.center
    if len(crosslink_bond_lengthList) == 1:
        r_max = 2 * crosslink_bond_lengthList[0] * (1 + bond_tolerance)
    else:
        r_max = sum(crosslink_bond_lengthList) * (1 + bond_tolerance)

    # Shift to freud's box frame once; all neighbor queries use this array.
    # TODO: Probably should shift all coordinates centered at 0 at beginning of
    # crosslink function
    shifted_coords = (coords - box_center).astype(np.float32)
    backbone_set = set(backbone_sites.tolist())
    n_backbones = len(backbone_sites)

    for site in rng.choice(np.arange(n_backbones), size=n_backbones, replace=False):
        origin = backbone_sites[site]

        if bond_graph is not None:  # validate distance along backbone chains
            excluded = set(
                nx.single_source_shortest_path_length(
                    bond_graph, origin, cutoff=minimum_intrachain_depth
                ).keys()
            )
        else:
            excluded = {origin}

        # _find_neighbors dispatches to AABBQuery (all-periodic) or tiled
        # cKDTree (any aperiodic axis).
        neighbors = _find_neighbors(
            shifted_coords,
            box,
            origin,
            r_max=r_max,
            r_min=1e-4,  # Will account for r_min during single site placement
        )

        for candidate_site in neighbors:
            if candidate_site not in backbone_set or candidate_site in excluded:
                continue
            yield (origin, int(candidate_site))


def _verify_one_site_bond(
    candidate,
    crosslinker,
    coords,
    minimum_separation,
    crosslink_bond_lengthsList,
    box,
    tolerance=0.1,
    rng=None,
    n_position_samples=8,
    n_rotation_samples=12,
):
    """Check bond_lengths, overlaps, and connectivity of a single_site crosslink."""
    unique_sites = crosslinker.unique_connection_sites
    if len(unique_sites) != 1:
        return False
    if crosslinker.n_connections != 2:
        return False

    pbc = box.pbc
    box_lengths = box.box_lengths  # Assume Cubic
    box_center = box.center
    connection_node = unique_sites[0]
    target_bond_lengths = np.asarray(crosslink_bond_lengthsList, dtype=float)
    if target_bond_lengths.size == 1:  # set for all connection sites
        target_bond_lengths = np.repeat(target_bond_lengths, crosslinker.n_connections)
    p1, p2 = coords[list(candidate)]
    r1, r2 = target_bond_lengths

    candidate_positions = _solve_two_bond_positions(
        p1, p2, r1, r2, tolerance, pbc, box_lengths, n_samples=n_position_samples
    )
    if candidate_positions is None:
        return False

    if rng is None:
        rng = np.random.default_rng()

    local_coords = crosslinker.coordinates - crosslinker.coordinates[connection_node]
    other_bead_mask = np.ones(crosslinker.n_sites, dtype=bool)
    other_bead_mask[connection_node] = False

    mask = np.ones(len(coords), dtype=bool)
    mask[list(candidate)] = False
    coords_excluding_candidate = coords[mask]

    # For wrapping for aperiodic axes
    aperiodic_axes = [ax for ax in range(3) if not pbc[ax]]
    box_lo = box_center - box_lengths / 2.0  # shape (3,)
    box_hi = box_center + box_lengths / 2.0

    checked_positions = set()
    for position in candidate_positions:
        positionTuple = tuple(position)  # skip duplicate sites
        if positionTuple not in checked_positions:
            checked_positions.add(positionTuple)
        else:
            continue
        for _ in range(n_rotation_samples):
            rotation = Rotation.random(random_state=rng)
            transformed = rotation.apply(local_coords) + position
            connection_pos = transformed[connection_node]

            if not _connection_bond_ok(
                connection_pos, p1, p2, r1, r2, tolerance, pbc, box_lengths
            ):
                continue
            if any(pbc):
                shifted = (transformed - box_center).astype(np.float32)
                newly_transformed = box.wrap(shifted).astype(np.float32) + box_center
                nan_mask = ~np.isnan(newly_transformed)
                transformed[nan_mask] = newly_transformed[nan_mask]
                connection_pos = transformed[connection_node]
            if aperiodic_axes:
                out_of_bounds = False
                for ax in aperiodic_axes:
                    if np.any(transformed[:, ax] < box_lo[ax] - 1e-5) or np.any(
                        transformed[:, ax] > box_hi[ax] + 1e-5
                    ):
                        out_of_bounds = True
                        break
                if out_of_bounds:
                    continue
            if not check_path(
                coords_excluding_candidate,
                connection_pos,
                radius=minimum_separation,
                tolerance=0.0,
                pbc=pbc,
                box_lengths=box_lengths,
            ):
                continue
            other_beads = transformed[other_bead_mask]
            if other_beads.shape[0] > 0 and not _check_all_overlap(
                other_beads, coords, minimum_separation, pbc, box_lengths
            ):
                continue

            return transformed

    return False


def _solve_two_bond_positions(p1, p2, r1, r2, tolerance, pbc, box_lengths, n_samples=8):
    """Analytic sphere-intersection solver, minimum-image aware. Same
    logic as before, but the p1->p2 displacement (and everything
    downstream of it) uses the minimum-image vector rather than the raw
    difference, since p1 and p2 may in fact be closer across a periodic
    face than their unwrapped coordinates suggest.
    """
    diff = _minimum_image(p2 - p1, pbc, box_lengths)
    d = np.linalg.norm(diff)
    if d < 1e-12:
        return None

    r1_lo, r1_hi = r1 * (1 - tolerance), r1 * (1 + tolerance)
    r2_lo, r2_hi = r2 * (1 - tolerance), r2 * (1 + tolerance)

    min_gap = max(0.0, r1_lo - r2_hi, r2_lo - r1_hi) - 1e-6
    max_reach = r1_hi + r2_hi
    if d < min_gap or d > max_reach:
        return None

    r1_eff, r2_eff = r1, r2
    if d > r1_eff + r2_eff:
        deficit = d - (r1_eff + r2_eff)
        headroom1, headroom2 = r1_hi - r1_eff, r2_hi - r2_eff
        r1_eff += deficit * (headroom1 / (headroom1 + headroom2))
        r2_eff += deficit * (headroom2 / (headroom1 + headroom2))
    elif d < abs(r1_eff - r2_eff):
        excess = abs(r1_eff - r2_eff) - d
        if r1_eff > r2_eff:
            room_down, room_up = r1_eff - r1_lo, r2_hi - r2_eff
            r1_eff -= excess * (room_down / (room_down + room_up))
            r2_eff += excess * (room_up / (room_down + room_up))
        else:
            room_down, room_up = r2_eff - r2_lo, r1_hi - r1_eff
            r2_eff -= excess * (room_down / (room_down + room_up))
            r1_eff += excess * (room_up / (room_down + room_up))

    u_hat = diff / d
    a = (d**2 + r1_eff**2 - r2_eff**2) / (2 * d)
    h_sq = r1_eff**2 - a**2
    if h_sq < 0:
        return None
    h = np.sqrt(h_sq)
    # Measured from p1 along the minimum-image direction toward p2 --
    # NOT p1's raw unwrapped distance to p2.
    circle_center = p1 + a * u_hat

    arbitrary = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(arbitrary, u_hat)) > 0.9:
        arbitrary = np.array([0.0, 1.0, 0.0])
    v1 = arbitrary - np.dot(arbitrary, u_hat) * u_hat
    v1 /= np.linalg.norm(v1)
    v2 = np.cross(u_hat, v1)

    thetas = np.linspace(0, 2 * np.pi, n_samples, endpoint=False)
    return circle_center + h * (
        np.outer(np.cos(thetas), v1) + np.outer(np.sin(thetas), v2)
    )


def _connection_bond_ok(position, p1, p2, r1, r2, tolerance, pbc, box_lengths):
    d1 = np.linalg.norm(_minimum_image(position - p1, pbc, box_lengths))
    d2 = np.linalg.norm(_minimum_image(position - p2, pbc, box_lengths))
    return r1 * (1 - tolerance) <= d1 <= r1 * (1 + tolerance) and r2 * (
        1 - tolerance
    ) <= d2 <= r2 * (1 + tolerance)


def _check_all_overlap(points, existing_points, radius, pbc, box_lengths):
    for p in points:
        if not check_path(
            existing_points,
            p,
            radius=radius,
            tolerance=0.0,
            pbc=pbc,
            box_lengths=box_lengths,
        ):
            return False
    return True


# =============================================================================
# Two-Site Crosslinker
# =============================================================================


def two_site_crosslink():
    """Crosslinking where backbones bond to two different crosslinker sites."""
    pass


# =============================================================================
# Periodic and neighbor utilities
# =============================================================================


def _minimum_image(diff, pbc, box_lengths):
    """Apply minimum-image convention to displacement vector(s) along
    periodic axes. Same convention as `check_path`'s inline wrap, kept
    as a shared helper so every distance calc here is consistent with it.

    Parameters
    ----------
    diff : np.ndarray (..., 3)
        Displacement vector(s), e.g. point - reference.
    pbc : np.ndarray (3,) bool
    box_lengths : np.ndarray (3,) float

    Returns
    -------
    np.ndarray, same shape as diff
    """
    diff = np.array(diff, dtype=float, copy=True)
    for axis in range(3):
        if pbc[axis] and not np.isinf(box_lengths[axis]):
            length = box_lengths[axis]
            diff[..., axis] -= np.round(diff[..., axis] / length) * length
    return diff


def _find_neighbors(coords, box, origin_idx, r_max, r_min=0.0):
    """Return neighbor indices of `origin_idx` within [r_min, r_max].

    Uses scikit cKDTree for handling MIC

    Parameters
    ----------
    coords    : np.ndarray (N, 3) float32, freud convention [-L/2, L/2).
    box : mbuild.path.constraints.CuboidConstraint  (periodic flags already set).
    origin_idx : int
    r_max      : float
    r_min      : float, default 0.0

    Returns
    -------
    np.ndarray of int  –  indices into `coords`.
    """

    periodic = box.pbc
    Ls = box.box_lengths
    boxsize = np.where(periodic & np.isfinite(Ls) & (Ls > 0), Ls, 0.0)

    pts = np.asarray(coords, dtype=np.float64).copy()
    for ax in range(3):  # cKDTree needs periodic coords in [0, L)
        if boxsize[ax] > 0:
            pts[:, ax] = np.mod(pts[:, ax], boxsize[ax])
            pts[pts[:, ax] >= boxsize[ax], ax] = 0.0  # np.mod can return exactly L

    tree = cKDTree(pts, boxsize=boxsize if np.any(boxsize > 0) else None)
    raw = np.asarray(tree.query_ball_point(pts[origin_idx], r=r_max), dtype=np.intp)
    if raw.size == 0:
        return np.empty(0, dtype=np.intp)

    d = pts[raw] - pts[origin_idx]
    for ax in range(3):
        if boxsize[ax] > 0:
            d[:, ax] -= np.round(d[:, ax] / boxsize[ax]) * boxsize[ax]
    dists = np.linalg.norm(d, axis=1)

    keep = (raw != origin_idx) & (dists >= r_min)
    return raw[keep]
