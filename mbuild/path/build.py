"""Classes to generate intra-molecular paths and configurations."""

import logging
import math
import time
from itertools import combinations_with_replacement

import networkx as nx
import numpy as np
from scipy.interpolate import interp1d

from mbuild import Box, Compound
from mbuild.exceptions import PathConvergenceError
from mbuild.path.constraints import CuboidConstraint, CylinderConstraint
from mbuild.path.namers import BEAD_NAME_DTYPE, BeadNamer
from mbuild.path.path_utils import (
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
from mbuild.utils.io import import_

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
    bead_name : str or numpy.ndarray (N, 1), default '_A'
        The name assigned to each site. This is helpful when using
        multiple `Path` instances to build heterogeneous systems.
        If an array of bead names is passed, it should be the same length
        as ``coordinates`` or the number of nodes in ``bond_graph``
        The array will be cast to BEAD_NAME_DTYPE, so bead names
        should not exceed MAX_BEAD_NAME_LENGTH characters.

    """

    def __init__(self, coordinates=None, bond_graph=None, bead_name="_A"):
        # Passing in an array of coordinates, bond graph and array of bead names
        # TODO: Cast to fp32 here?
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
            self.beads = bead_name.astype(BEAD_NAME_DTYPE)
        # Passing in an array of coordinates, bond graph, and single bead name
        elif coordinates is not None and bond_graph is not None:
            assert len(coordinates) == len(bond_graph)
            self.bond_graph = bond_graph
            self.coordinates = coordinates
            self.beads = np.array(
                [str(bead_name) for _ in range(len(coordinates))], dtype=BEAD_NAME_DTYPE
            )
        # Only passing in a bond graph with data defined for xyz and name
        # TODO: Cast to fp32 here?
        elif coordinates is None and bond_graph is not None:
            self.bond_graph = bond_graph
            self.coordinates = np.array(  # might not have a bond_graph xyz
                [node.get("xyz") for node in bond_graph.nodes(data=True)]
            )
            self.beads = np.array(
                [node.get("name") for node in bond_graph.nodes(data=True)]
            ).astype(BEAD_NAME_DTYPE)
        # Only passing in coordinates, need to create bond graph object.
        elif coordinates is not None and bond_graph is None:
            self.coordinates = np.asarray(coordinates)
            self.bond_graph = nx.Graph()
            # Only passed a single string for bead name, create array
            if isinstance(bead_name, str):
                self.beads = np.array(
                    [bead_name for _ in coordinates], dtype=BEAD_NAME_DTYPE
                )
            # Passed array of bead names, cast to the bead name dtype
            elif isinstance(bead_name, np.ndarray):
                self.beads = bead_name.astype(BEAD_NAME_DTYPE)
            elif isinstance(bead_name, list):
                self.beads = np.array(bead_name, dtype=BEAD_NAME_DTYPE)
            for idx in range(len((self.coordinates))):
                self.bond_graph.add_node(idx)
        # Nothing is defined, create empty place holders for coords, bond graph and bead names
        else:
            self.coordinates = np.array([], dtype=np.float32)
            self.bond_graph = nx.Graph()
            self.beads = np.array([], dtype=BEAD_NAME_DTYPE)

    def __eq__(self, other):
        return (
            np.all(self.coordinates == other.coordinates)
            and np.all(self.beads == other.beads)
            and nx.is_isomorphic(self.bond_graph, other.bond_graph)
        )

    def __add__(self, other):
        coordinates = np.concat((self.coordinates, other.coordinates))
        beads = np.concat((self.beads, other.beads))
        offset = len(self)
        mapping = {node: node + offset for node in range(len(other))}
        shifted_graph = nx.relabel_nodes(other.bond_graph, mapping)
        bond_graph = nx.compose(self.bond_graph, shifted_graph)
        return Path(coordinates, bond_graph, beads)

    def __len__(self):
        if hasattr(self, "coordinates"):
            return len(self.coordinates)
        return 0

    def add_path(self, otherPath, other_coordinates=None):
        """Add one path to current path."""
        if other_coordinates is None:
            other_coordinates = otherPath.coordinates
        initial_n_nodes = len(self)
        self.append_coordinates(other_coordinates, otherPath.beads)
        for site1, site2 in otherPath.bond_graph.edges():
            self.bond_graph.add_edge(site1 + initial_n_nodes, site2 + initial_n_nodes)

    @classmethod
    def from_compound(cls, compound):
        coordinates = compound.xyz
        names = [particle.name for particle in compound.particles()]

        # Create the path with coordinates and bond graph
        path = cls(coordinates=coordinates, bead_name=names)
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

    def append_coordinates(self, points, bead_names):
        """Create new coordinates for appending values to.

        Parameters
        ----------
        points : array-like of float, shape (N, 3), required
            A set of 3D coordinates to append to Path.coordinates.
        names : array-like of str, shape (N, 3), required
            The set of bead names that correspond to points, appends to Path.beads.
        """
        if hasattr(points, "__len__") and len(points) == 3 and points.ndim == 1:
            points = np.asarray([points])
        if not isinstance(points, np.ndarray):
            try:
                points = np.asarray(points)
            except TypeError:
                raise ValueError(f"{points=} must be an array of shape (N,3)")
        if points.ndim == 1:
            points = np.array([points])  # make a 2d array
        # Create sequence of bead names
        if isinstance(bead_names, str):
            bead_names = np.array([bead_names] * len(points), dtype=BEAD_NAME_DTYPE)

        if self.coordinates.size == 0:
            self.coordinates = points
            self.beads = bead_names
            self._extend_bond_graph()
            return
        # Append to previous set of coordinates and names
        self.coordinates = np.concatenate((self.coordinates, points))
        self.beads = np.concatenate((self.beads, bead_names))
        self._extend_bond_graph()

    def _extend_coordinates(self, N):
        """Create new coordinate and bead name place holders for setting values."""
        if self.coordinates.size == 0:
            self.coordinates = np.zeros((N, 3), dtype=np.float32)
            self.beads = np.zeros(N, dtype=BEAD_NAME_DTYPE)  # Place holder is empty str
            return
        # Update coordinates array
        zeros = np.zeros((N, 3), dtype=self.coordinates.dtype)
        self.coordinates = np.concatenate([self.coordinates, zeros])
        # Update bead names array
        empty = np.zeros(N, dtype=BEAD_NAME_DTYPE)
        self.beads = np.concatenate([self.beads, empty])

    def _extend_bond_graph(self):
        """Make sure all coordinates are properly added to the bondgraph."""
        if len(self.bond_graph) == len(self.coordinates):
            return
        index = len(self.bond_graph)
        for _ in self.coordinates[len(self.bond_graph) :]:
            name = self.beads[index]
            self.bond_graph.add_node(index, name=name)
            index += 1

    def _extend_beads(self, bead_name):
        """Update bead names to extended coordinates."""
        diff = len(self.coordinates) - len(self.beads)
        if not diff:
            return
        self.beads = np.concatenate((self.beads, [bead_name] * diff))

    # TODO: Can we accept negative attach indices?
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
        self.bond_graph.add_edge(
            u_of_edge=u,
            v_of_edge=v,
            direction=bond_vec.tolist(),
            length=float(bond_length),
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
        # TODO: Should we have a mass parameter? Could be useful for density termination
        # TODO: Could also add an is_atomistic flag and validate the bead names here before sending to Compound.
        for node_id in self.bond_graph.nodes:
            compounds.append(
                Compound(
                    name=self.beads[node_id], pos=self.coordinates[node_id], mass=1.0
                )
            )
        compound.add(compounds)
        for edge1, edge2 in self.bond_graph.edges():
            compound.bond_graph.add_edge(
                compounds[edge1], compounds[edge2], bond_order=1.0
            )
        return compound

    def to_cgsmiles_graph(self, fragname_map=None):
        """Convert this path's bond graph to a CGsmiles-compatible meta graph.

        Parameters
        ----------
        fragname_map : dict[str, str], optional
            Mapping of bead names to CGsmiles fragment names. Bead names not
            present in the map are used as fragment names directly. Useful
            when path bead names differ from the fragment
            names used in the CGsmiles fragment string.

        Returns
        -------
        networkx.Graph
            Meta graph usable with ``cgsmiles.MoleculeResolver``.

        See ``mbuild.coarse_graining.to_cgsmiles_graph``.
        """
        from mbuild.coarse_graining import to_cgsmiles_graph

        return to_cgsmiles_graph(self, fragname_map=fragname_map)

    def to_cgsmiles(self, fragname_map=None):
        """Write the coarse-grained level of this Path as a CGsmiles string.

        The returned string describes the bead sequence and connectivity
        (including branches and rings) at the coarse-grained level.
        Append fragment definitions (e.g. ``"{#A=[>]CC[<]}"``)
        to obtain a fully resolvable CGsmiles string.

        Parameters
        ----------
        fragname_map : dict[str, str], optional
            Mapping of bead names to CGsmiles fragment names. Bead names not
            present in the map are used as fragment names directly. Useful
            when path bead names differ from the fragment
            names used in the CGsmiles fragment string.

        Returns
        -------
        str
            The CGsmiles graph string, e.g. ``"{[#A][#A]([#B][#B])[#A]}"``.
            Paths holding multiple disconnected molecules (e.g. a box of
            chains) are written as ``.``-separated segments, which CGsmiles
            reads as zero-order (non-bonded) connections.

        See ``mbuild.coarse_graining.to_cgsmiles``.
        """
        from mbuild.coarse_graining import to_cgsmiles

        return to_cgsmiles(self, fragname_map=fragname_map)

    def backmap(self, fragments=None, **kwargs):
        """Backmap the path to an atomistic Compound using CGsmiles.

        Resolves each bead to molecular detail and returns an atomistic
        ``mbuild.Compound`` that retains the path's conformation. Works
        for any bond graph topology, including branch points. Fragments
        are defined by CGsmiles fragment strings (SMILES with bonding
        descriptors), by tagged mBuild compounds passed via
        ``templates``, or a mix of both.

        See ``mbuild.coarse_graining.backmap`` for parameters.

        Example
        -------
        >>> path = straight_line(spacing=0.25, N=10, bead_name="PEO")
        >>> compound = path.backmap("{#PEO=[>]COC[<]}")
        >>> # or, defining the fragment with a tagged compound instead
        >>> template = mb.load("C{>}O{ }C{<}", smiles=True)  # doctest: +SKIP
        >>> compound = path.backmap(templates={"PEO": template})  # doctest: +SKIP
        """
        from mbuild.coarse_graining import backmap

        return backmap(self, fragments, **kwargs)

    def to_mol2(self):
        """Convert a path to a .mol2 file."""
        from mbuild.path.formats import to_mol2 as _to_mol2

        return _to_mol2(self)

    def to_mol(self):
        """Convert mBuild Path to SDF/MOL format."""
        from mbuild.path.formats import to_mol as _to_mol

        return _to_mol(self)

    def to_mol3000(self):
        """Convert mBuild Path to SDF/MOL V3000 format."""
        from mbuild.path.formats import to_mol3000 as _to_mol3000

        return _to_mol3000(self)

    def visualize(self, radius, hide_periodic_bonds=False):
        """Visualize a path with py3Dmol.

        Parameters
        ----------
        radius : float, required
            Sets the bead size in py3Dmol.
        hide_periodic_bonds : bool, default False
            If ``True`` bonds crossing periodic boundaries aren't shown.
        """
        from mbuild.utils.visualize import visualize_path

        return visualize_path(self, radius, hide_periodic_bonds)

    def relax(
        self,
        bead_radius,
        bond_length=None,
        angles_sampler=None,
        steps=1000,
        seed=1,
        run_on_gpu=False,
    ):
        """Relax the path with a coarse Kremer-Grest forcefield.

        Builds a WCA + FENE forcefield from the given parameters and runs a
        capped-displacement warm-up followed by FIRE minimization in HOOMD.
        Relaxed coordinates are synced back onto the path.

        Parameters
        ----------
        bead_radius : float or dict
            WCA bead radius (center-to-center exclusion diameter). A dict keyed
            by bead type gives each bead type its own size.
        bond_length : float or dict, optional
            Target bond length. Defaults to per-bond-type means over the path.
            A dict is keyed by an unordered bead-type pair.
        angles_sampler : mbuild.path.points.AnglesSampler or dict, optional
            If given, adds a tabulated angle potential. A dict keyed by a
            bead-type triple parameterizes angles per type.
        steps : int, optional, default 1,000
            Number of FIRE minimization steps.
        seed : int, optional, default 1
            Random seed for the simulation.
        run_on_gpu : bool, default False
            If True, HOOMD will attempt to run on the GPU if the user has
            a GPU compatible version of HOOMD installed, and a compatible
            GPU available. If one isn't found, it will revert to the CPU.
            if False, HOOMD will run on the CPU.
        """
        import mbuild.simulation
        from mbuild.utils.simulation.path_forces import (
            PathForcefield,
            _mean_bond_length,
        )

        # Scalar length scale for the displacement cap, even if per-type lengths
        # were given.
        if isinstance(bond_length, (int, float)):
            bond_eff = bond_length
        else:
            bond_eff = _mean_bond_length(self)
        forcefield = PathForcefield(
            radius=bead_radius, bond_length=bond_length, angles=angles_sampler
        )
        sim = mbuild.simulation.HoomdSimulation(
            self, forcefield=forcefield, seed=seed, run_on_gpu=run_on_gpu
        )
        sim.cap_displacement(n_steps=500, dt=1, max_displacement=0.01 * bond_eff)
        sim.fire(n_steps=steps)

    def print_bond_lengths(self, box=None):
        """Compute and return info about path bonds.

        Parameters
        ----------
        box : list, optional
            list of [Lx, Ly, Lz] to use for checking for periodic bonds. Assume centered at (0,0,0)

        Returns
        -------
        bonds : dict
            Dictionary mapping (i, center, k) tuples to their angles in nm.
        bonds_typesDict : dict
            Dictionary mapping sorted bead-type tuples to lists of bonds.
        """
        positions = self.coordinates
        bond_lengths = {}
        beads = set(self.beads)
        bond_types = list(combinations_with_replacement(beads, 2))
        bond_typesDict = {
            tuple(sorted((str(b1), str(b2)))): [] for b1, b2 in bond_types
        }

        for i, j in self.bond_graph.edges():
            delta = positions[j] - positions[i]
            if isinstance(box, CuboidConstraint):
                box_arr = box.box_lengths
                delta -= np.round(delta / box_arr) * box_arr
            elif isinstance(box, Box):
                box_arr = box.lengths
                delta -= np.round(delta / box_arr) * box_arr
            elif isinstance(box, list):
                box_arr = np.array(box)
                delta -= np.round(delta / box_arr) * box_arr
            bl = np.linalg.norm(delta)
            bond_lengths[(i, j)] = bl
            bond_name = tuple(sorted((str(self.beads[i]), str(self.beads[j]))))
            bond_typesDict[bond_name].append(bl)

        # Summary stats
        print(f"Min bond length: {min(bond_lengths.values()):.4f}")
        print(f"Max bond length: {max(bond_lengths.values()):.4f}")
        print(f"Mean bond length: {np.mean(list(bond_lengths.values())):.4f}")
        for key in bond_typesDict:
            if not len(bond_typesDict[key]):
                continue
            print(
                f"{key}: Max={max(bond_typesDict[key]):.2f}, Min={min(bond_typesDict[key]):.2f}"
            )

        return bond_lengths, bond_typesDict

    def print_angle_lengths(self, box=None):
        """Compute and return info about path angles.

        Parameters
        ----------
        box : array-like of shape (3,), optional
            Box dimensions [Lx, Ly, Lz]. If None, no periodic wrapping is applied.

        Returns
        -------
        angles : dict
            Dictionary mapping (i, center, k) tuples to their angles in degrees.
        angle_typesDict : dict
            Dictionary mapping sorted bead-type tuples to lists of angles.
        """
        from itertools import combinations, combinations_with_replacement

        positions = self.coordinates
        angles = {}
        beads = set(self.beads)

        # Build angle type keys: (outer, center, outer) with outer pair sorted
        angle_types = set()
        for center_bead in beads:
            for b1, b2 in combinations_with_replacement(beads, 2):
                angle_types.add((str(b1), str(center_bead), str(b2)))
        angle_typesDict = {key: [] for key in angle_types}

        for center_node in self.bond_graph.nodes():
            neighbors = list(self.bond_graph.neighbors(center_node))
            if len(neighbors) < 2:
                continue

            for i_node, k_node in combinations(neighbors, 2):
                # Vector from center to i
                delta_i = positions[i_node] - positions[center_node]
                if box is not None:
                    box_arr = np.array(box)
                    delta_i -= np.round(delta_i / box_arr) * box_arr

                # Vector from center to k
                delta_k = positions[k_node] - positions[center_node]
                if box is not None:
                    delta_k -= np.round(delta_k / box_arr) * box_arr

                # Compute angle via dot product
                cos_angle = np.dot(delta_i, delta_k) / (
                    np.linalg.norm(delta_i) * np.linalg.norm(delta_k)
                )
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angle_deg = np.degrees(np.arccos(cos_angle))

                angles[(i_node, center_node, k_node)] = angle_deg

                # Categorize by type (sort outer beads for consistency)
                outer_beads = tuple(
                    sorted((str(self.beads[i_node]), str(self.beads[k_node])))
                )
                center_bead = str(self.beads[center_node])
                angle_key = (outer_beads[0], center_bead, outer_beads[1])
                angle_typesDict[angle_key].append(angle_deg)

        # Summary stats
        print(f"Min angle: {min(angles.values()):.2f}°")
        print(f"Max angle: {max(angles.values()):.2f}°")
        print(f"Mean angle: {np.mean(list(angles.values())):.2f}°")
        for key in angle_typesDict:
            if not len(angle_typesDict[key]):
                continue
            print(
                f"{key}: Max={max(angle_typesDict[key]):.2f}°, Min={min(angle_typesDict[key]):.2f}°, Mean={np.mean(angle_typesDict[key]):.2f}°"
            )

        return angles, angle_typesDict

    def freud_rdf(self, box=None, bins=50, r_max=1):
        freud = import_("freud")
        if box is None:
            box = np.array([np.inf, np.inf, np.inf, 0, 0, 0])  # infinite box
        elif isinstance(box, Box):
            box = np.array([*box.lengths, 0, 0, 0])
        elif len(box) == 3:
            box = np.array([*box, 0, 0, 0])  # assume orthorhombic

        rdf = freud.density.RDF(bins=bins, r_max=r_max)
        rdf.compute(system=(box, self.coordinates))
        return rdf


def lamellar(
    path=None,
    num_layers=1,
    layer_separation=None,
    layer_length=None,
    spacing=None,
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
    spacing : float (nm), required
        The distance between two adjacent sites in the path.
    initial_point : nd.array (1,3), default (0,0,0)
        The coordinate of the first site of the lamellar path.
    num_stacks : int, default 1
        The number of times to repeat each layer in the Z direction.
    stack_separation : float (nm), optional
        The distance between two stacked layers. Required if `num_stacks` >= 2.
    left_to_right : boolean, default True
        If `True`, the first layer is built with increasing y-coordinates from the origin.
    bead_name : str or BeadNamer, optional, default '_A'
        Name(s) to assign to beads. A plain string assigns the same name to
        every bead. Pass a ``BeadNamer`` instance (e.g. ``CyclicNamer``,
        ``RandomNamer``, ``MarkovNamer``) for heterogeneous sequences.
    """
    if path is None:
        path = Path()
    initial_point = np.asarray(initial_point)

    # Coordinates in the y-direction (layer-length) of the lamellar layer
    layer_spacing = np.arange(0, layer_length, spacing)
    if not left_to_right:
        layer_spacing *= -1
    layer_spacing += initial_point[1]

    # Info needed for generating coords of the arc curves between layers
    r = layer_separation / 2
    arc_length = r * np.pi
    arc_num_points = math.floor(arc_length / spacing)
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
        arc_num_points = math.floor(arc_length / spacing)
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
    # Coordinates are set, update the Path object
    start_index = len(path.coordinates)
    stop_index = start_index + len(coordinates)
    namer = BeadNamer.coerce(bead_name)
    names = np.array(
        [next(namer) for _ in range(len(coordinates))], dtype=BEAD_NAME_DTYPE
    )
    path.append_coordinates(coordinates, names)
    path._connect_edges(
        connectivity="linear", indices=np.arange(start_index, stop_index)
    )
    return path


def straight_line(spacing, N, path=None, direction=(1, 0, 0), bead_name="_A"):
    """Generates a set of coordinates in a straight line along a given axis.

    Parameters
    ----------
    N : int, required
        The number of sites in the path.
    spacing : float, required
        The distance between sites along the path.
    path : mbuild.path.Path
        The Path object to populate with coordinates.
        Defaults to creating a new, empty Path if no Path is given.
    direction : array-like (1,3), default = (1,0,0)
        The direction to align the straight path along.
    bead_name : str or path.BeadNamer, optional, default '_A'
        Name(s) to assign to beads. A plain string assigns the same name to
        every bead. Pass a ``BeadNamer`` instance for heterogeneous sequences.
        See mbuild.path.namers.py
    """
    if path is None:
        path = Path()
    direction = np.asarray(direction)
    coordinates = np.array([np.zeros(3) + i * spacing * direction for i in range(N)])
    start_index = len(path.coordinates)
    stop_index = start_index + N
    namer = BeadNamer.coerce(bead_name)
    names = np.array(
        [next(namer) for _ in range(len(coordinates))], dtype=BEAD_NAME_DTYPE
    )
    path.append_coordinates(coordinates, names)
    path._connect_edges(
        connectivity="linear", indices=np.arange(start_index, stop_index)
    )
    return path


def cyclic(spacing=None, N=None, path=None, radius=None, closed=True, bead_name="_A"):
    """Generates a set of coordinates evenly spaced along a circle.

    Parameters
    ----------
    spacing : float, optional
        Distance between sites along the path.
    N : int, optional
        Number of sites in the cyclic path.
    path : mbuild.path.Path
        The Path object to populate with coordinates.
        Defaults to creating a new, empty Path if no Path is given.
    radius : float, optional
        The radius (nm) of the cyclic path.
    closed : bool, default True
        If `True` the cyclic path is closed by bonding the first and last sites together
    bead_name : str or path.BeadNamer, optional, default '_A'
        Name(s) to assign to beads. A plain string assigns the same name to
        every bead. Pass a ``BeadNamer`` instance for heterogeneous sequences.
        See mbuild.path.namers.py

    Notes
    -----
    Only two of spacing, N and radius can be defined, as the third
    is determined by the other two.
    """
    if path is None:
        path = Path()
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
    coordinates = np.array(
        [(np.cos(a) * radius, np.sin(a) * radius, 0) for a in angles]
    )
    start_index = len(path.coordinates)
    stop_index = start_index + len(coordinates)
    namer = BeadNamer.coerce(bead_name)
    names = np.array(
        [next(namer) for _ in range(len(coordinates))], dtype=BEAD_NAME_DTYPE
    )
    path.append_coordinates(coordinates, names)
    if closed:
        path._connect_edges(
            connectivity="cycle", indices=np.arange(start_index, stop_index)
        )
    else:
        path._connect_edges(
            connectivity="linear", indices=np.arange(start_index, stop_index)
        )
    return path


def knot(spacing, N, m, path=None, closed=True, bead_name="_A"):
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
    path : mbuild.path.Path
        The Path object to populate with coordinates.
        Defaults to creating a new, empty Path if no Path is given.
    closed : bool, default True
        If `True` the cyclic path is closed by bonding the first and last sites together
    bead_name : str or path.BeadNamer, optional, default '_A'
        Name(s) to assign to beads. A plain string assigns the same name to
        every bead. Pass a ``BeadNamer`` instance for heterogeneous sequences.
        See mbuild.path.namers.py
    """
    if path is None:
        path = Path()
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
    coordinates = np.stack((x_interp, y_interp, z_interp), axis=1)

    start_index = len(path.coordinates)
    stop_index = start_index + len(coordinates)
    namer = BeadNamer.coerce(bead_name)
    names = np.array(
        [next(namer) for _ in range(len(coordinates))], dtype=BEAD_NAME_DTYPE
    )
    path.append_coordinates(coordinates, names)
    if closed:
        path._connect_edges(
            connectivity="cycle", indices=np.arange(start_index, stop_index)
        )
    else:
        path._connect_edges(
            connectivit="linear", indices=np.arange(start_index, stop_index)
        )
    return path


def helix(
    N, radius, rise, twist, path=None, right_handed=True, bottom_up=True, bead_name="_A"
):
    """Generate helical path.

    Parameters:
    -----------
    N : int, required
        Number of sites in the path
    radius : float, required
        Radius of the helix (nm)
    rise : float, required
        Rise per site on path (nm)
    twist : float, required
        Twist per site in path (degrees)
    path : mbuild.path.Path
        The Path object to populate with coordinates.
        Defaults to creating a new, empty Path if no Path is given.
    right_handed : bool, default True
        Set the handedness of the helical twist
    bottom_up : bool, default True
        If True, the twist is in the positive Z direction
    bead_name : str or path.BeadNamer, optional, default '_A'
        Name(s) to assign to beads. A plain string assigns the same name to
        every bead. Pass a ``BeadNamer`` instance for heterogeneous sequences.
        See mbuild.path.namers.py
    """
    if path is None:
        path = Path()
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

    start_index = len(path.coordinates)
    stop_index = start_index + len(coordinates)
    namer = BeadNamer.coerce(bead_name)
    names = np.array(
        [next(namer) for _ in range(len(coordinates))], dtype=BEAD_NAME_DTYPE
    )
    path.append_coordinates(coordinates, names)
    path._connect_edges(
        connectivity="linear", indices=np.arange(start_index, stop_index)
    )
    return path


def spiral_2D(N, a, b, spacing, path=None, bead_name="_A"):
    """Generate a 2D spiral path in the XY plane.

    Parameters
    ----------
    N : int, required
        Number of sites in the path
    a : float, required
        The initial radius (nm)
    b : float, required
        Determines how fast radius grows per angle increment (r = a + bθ)
    spacing : float, required
        Distance between adjacent sites (nm)
    path : mbuild.path.Path
        The Path object to populate with coordinates.
        Defaults to creating a new, empty Path if no Path is given.
    bead_name : str or path.BeadNamer, optional, default '_A'
        Name(s) to assign to beads. A plain string assigns the same name to
        every bead. Pass a ``BeadNamer`` instance for heterogeneous sequences.
        See mbuild.path.namers.py
    """
    if path is None:
        path = Path()
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

    start_index = len(path.coordinates)
    stop_index = start_index + len(coordinates)
    namer = BeadNamer.coerce(bead_name)
    names = np.array(
        [next(namer) for _ in range(len(coordinates))], dtype=BEAD_NAME_DTYPE
    )
    path.append_coordinates(coordinates, names)
    path._connect_edges(
        connectivity="linear", indices=np.arange(start_index, stop_index)
    )
    return path


def zigzag(
    N,
    spacing,
    path=None,
    angle_deg=120.0,
    sites_per_segment=4,
    plane="xy",
    bead_name="_A",
):
    """Generates a path following a zig-zag pattern in a given plane.

    Parameters
    ----------
    N : int, required
        Number of sites in the path
    spacing : float, required
        The distance (nm) between consecutive sites along the path.
    path : mbuild.path.Path
        The Path object to populate with coordinates.
        Defaults to creating a new, empty Path if no Path is given.
    angle_deg : float, default = 120.
        The rotation (degrees) applied between segments
    sites_per_segment : int, default = 4
        The number of sites before rotating and beginning next segment.
    plane : str, default = "xy"
        The plane that the sites in the path occupy
    bead_name : str or path.BeadNamer, optional, default '_A'
        Name(s) to assign to beads. A plain string assigns the same name to
        every bead. Pass a ``BeadNamer`` instance for heterogeneous sequences.
        See mbuild.path.namers.py
    """
    if N % sites_per_segment != 0:
        raise ValueError("N must be evenly divisible by sites_per_segment")

    if path is None:
        path = Path()
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

    start_index = len(path.coordinates)
    stop_index = start_index + len(coordinates)
    namer = BeadNamer.coerce(bead_name)
    names = np.array(
        [next(namer) for _ in range(len(coordinates))], dtype=BEAD_NAME_DTYPE
    )
    path.append_coordinates(coordinates, names)
    path._connect_edges(
        connectivity="linear", indices=np.arange(start_index, stop_index)
    )
    return path


def hard_sphere_random_walk(
    path=None,
    bead_name="_A",
    bond_length=0.15,
    radius=0.1,
    rw_angles=None,
    termination=None,
    volume_constraint=None,
    bias=None,
    include_compound=None,
    connectivity="linear",
    initial_point=None,
    seed=42,
    trial_batch_size=20,
    tolerance=1e-5,
    chunk_size=512,
):
    """Generates coordinates from a self avoiding random walk using
    fixed bond lengths, hard spheres, and minimum and maximum angles
    formed by 3 consecutive points.

    Parameters:
    -----------
    path : mbuild.path.Path
        The Path object to populate with coordinates.
        Defaults to creating a new, empty Path if no Path is given.
    bead_name : str or path.BeadNamer, optional, default '_A'
        Name(s) to assign to beads. A plain string assigns the same name to
        every bead. Pass a ``BeadNamer`` instance for heterogeneous sequences.
        See mbuild.path.namers.py
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
    include_compound : mbuild.compound.Compound, default None
        If an mBuild Compound is given, the random walk with include its coordinates
        when checking for overlapping sites.
    connectivity : str, default "linear"
        Will be used to connect the bond_graph of the walk post generation based on a specified method.
        See path._connect_edges for different options.
    initial_point : array-like or int, optional
        Used as the coordinate for the first site in this random walk path. If an integer is
        passed, look in coordinates of passed path object, and grab the starting coordinates from there.
    seed : int, default = 42
        Random seed
    trial_batch_size : int, default = 20
        The number of trial moves to attempt in parallel for each step.
    tolerance : float, default = 1e-5
        Tolerance used for rounding and checking for overlaps.
    chunk_size : int, default = 512
        Size of coordinate chunks to allocate
    """
    # Create seed sequence used by multiple path classes
    # The namer seed is separate, so that coordinates are impacted by naming methods.
    previous_count = len(path.coordinates) if path else 0
    seed_sequence = np.random.SeedSequence(seed + previous_count)
    name_seed_sequence = seed_sequence.spawn(1)[0]
    rng = np.random.default_rng(seed_sequence)

    # Create state object to track random walk progress
    state = RandomWalkState(
        bond_length=bond_length,
        radius=radius,
        angles_sampler=rw_angles,
        bead_name=bead_name,
        initial_point=initial_point,
        previous_count=previous_count,
        include_compound=include_compound,
        connectivity=connectivity,
        seed=seed,
        volume_constraint=volume_constraint,
        tolerance=tolerance,
        trial_batch_size=int(trial_batch_size),
        chunk_size=chunk_size,
        rng=rng,
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

    namer = BeadNamer.coerce(bead_name)
    namer._attach_rng(np.random.default_rng(name_seed_sequence))

    # Set up PBC info from volume constraints
    # TODO: We can probably out-source pbc, box_lengths return to the Constraint classes
    if isinstance(volume_constraint, CuboidConstraint):
        pbc = np.asarray(volume_constraint.pbc, dtype=np.bool_)
        box_lengths = volume_constraint.box_lengths.astype(np.float32)
    elif isinstance(volume_constraint, CylinderConstraint):
        pbc = np.array(
            [False, False, volume_constraint.periodic_height], dtype=np.bool_
        )
        box_lengths = np.array(
            [
                volume_constraint.radius * 2,
                volume_constraint.radius * 2,
                volume_constraint.height,
            ]
        ).astype(np.float32)
    else:
        pbc = np.array([False, False, False], dtype=np.bool_)
        box_lengths = np.array([np.inf, np.inf, np.inf], dtype=np.float32)
    state.pbc = pbc
    state.box_lengths = box_lengths

    # Set up bias conditions
    if bias:
        bias._attach_path(path, state)
        state.bias = bias

    # Initialize coordinates based on starting conditions
    # The path used for this random walk contains previous sites
    if not path.coordinates.size == 0:
        coordinates = np.concatenate(
            (
                path.coordinates.astype(np.float32),
                np.zeros((chunk_size, 3), dtype=np.float32),
            ),
            axis=0,
        )
        beads = np.concatenate(
            (path.beads, np.zeros(chunk_size, dtype=BEAD_NAME_DTYPE)),
            axis=0,
        )
        state.count = len(path.coordinates)  # starting index
    # The path used for this RW doesn't have previous sites
    else:
        coordinates = np.zeros((chunk_size, 3), dtype=np.float32)
        beads = np.zeros(chunk_size, dtype=BEAD_NAME_DTYPE)
        state.count = 0

    state.init_count = state.count

    # Select methods for random walk
    check_path_cpu = check_path
    next_step = random_coordinate

    # Set start time for wall time terminator
    state.start_time = time.time()

    # Set the first two coordinates
    for num_tries in range(100):
        # Get the first coordinate
        initial_xyz = get_initial_point(
            state=state,
            existing_points=coordinates[: state.count],
            beads=beads[: state.count],
            check_path=check_path_cpu,
            next_step=next_step,
        )
        coordinates[state.count] = initial_xyz
        beads[state.count] = next(namer)
        state.count += 1
        if state.check_termination(path, coordinates, beads):
            return path

        # Get the second coordinate
        second_xyz = get_second_point(
            state=state,
            existing_points=coordinates[: state.count],
            beads=beads[: state.count],
            check_path=check_path_cpu,
            next_step=next_step,
        )
        if second_xyz is not None:
            coordinates[state.count] = second_xyz
            beads[state.count] = next(namer)
            state.count += 1
            break
        else:
            # Retrying first and second point, reduce count by 1
            state.count -= 1
            if num_tries == 99:
                if state.initial_point is not None:
                    raise PathConvergenceError(
                        f"Failed to initiate random walk with {initial_point=}. Try a different initial_point."
                    )
                else:
                    raise PathConvergenceError(
                        f"Failed after {num_tries + 1} to generate a starting point. System is probably too densely packed."
                    )

    if state.check_termination(path, coordinates, beads):
        return path

    # Main random walk loop
    walk_finished = False
    while not walk_finished:
        batch_angles, batch_vectors = generate_trials(state)
        candidates = next_step(
            pos1=coordinates[state.count - 1],
            pos2=coordinates[state.count - 2],
            bond_length=bond_length,
            thetas=batch_angles,
            r_vectors=batch_vectors,
        )
        if state.volume_constraint:
            is_inside_mask = volume_constraint.is_inside(
                points=candidates, buffer=tolerance
            )
            candidates = candidates[is_inside_mask]
        # If there is a bias, sort candidates according to the bias
        if state.bias:
            candidates = bias(
                candidates=candidates,
                coordinates=coordinates[: state.count],
                names=beads[: state.count],
            )
        # Handle postion for PBCs.
        if any(pbc):
            candidates = (
                volume_constraint.mins
                + np.mod(candidates - volume_constraint.mins, box_lengths)
            ).astype(np.float32)
        # Check candidate sites
        accept_xyz = None
        existing_points = coordinates[: state.count]
        if state.include_compound:  # Include compound's particle coordinates
            existing_points = np.concat((existing_points, include_compound.xyz))
        # Iterate through current state of candidates, break after first accept
        for xyz in candidates:
            if check_path_cpu(
                existing_points=existing_points,
                new_point=xyz,
                radius=radius,
                tolerance=tolerance,
                pbc=pbc,
                box_lengths=box_lengths,
            ):
                accept_xyz = xyz
                break

        if accept_xyz is not None:
            coordinates[state.count] = accept_xyz
            beads[state.count] = next(namer)
            state.count += 1
        state.attempts += 1

        # Extend coordinates array if we're running out of space
        if state.count + 1 >= len(coordinates):
            path.coordinates = coordinates  # Save progress first
            path.beads = beads
            path._extend_coordinates(N=chunk_size)
            coordinates = path.coordinates
            beads = path.beads

        walk_finished = termination.is_met(
            coordinates=coordinates[: state.count], names=beads[: state.count]
        )
    state.check_termination(path, coordinates, beads)

    return path


def _normalize_initial_point(initial_point):
    """Resolve a user given initial_point into a value and what it means.

    An initial point is either a coordinate to start the walk at, or an index
    of a site in an existing path to start the walk from. A sequence of 3
    values is a coordinate. Any single integer is an index, including numpy
    integer scalars and single element integer arrays such as ``np.int64(4)``,
    ``np.array(4)`` and ``np.array([4])``. Indices are returned as python ints.

    Parameters
    ----------
    initial_point : array-like, int or None
        The initial_point argument given to a random walk.

    Returns
    -------
    value : np.ndarray of shape (3,), int or None
        The coordinate, the site index, or None if no initial point was given.
    starting_from_site : bool
        True when value is an index of a site in an existing path.

    Raises
    ------
    ValueError
        If initial_point is neither a 3 coordinate array nor an integer index.
    """
    if initial_point is None:
        return None, False
    point = np.asarray(initial_point)
    if point.ndim == 1 and point.size == 3:
        return point, False
    if point.size == 1 and np.issubdtype(point.dtype, np.integer):
        return point.item(), True
    raise ValueError(
        f"Unsupported initial_point {initial_point!r}. Pass either an "
        "array-like of 3 coordinates, or an integer index of a site in the "
        "coordinates of an existing path."
    )


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
    initial_point : np.ndarray, int or None
        Specified initial coordinate, or the index of the site in an existing
        path that this walk starts from
    starting_from_site : bool
        True when initial_point is an index of a site in an existing path
    include_compound : mbuild.compound.Compound, default None
        If an mBuild Compound is given, the random walk with include its coordinates
        when checking for overlapping sites.
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
        include_compound=None,
        seed=42,
        volume_constraint=None,
        termination=None,
        bias=None,
        tolerance=1e-5,
        trial_batch_size=20,
        chunk_size=512,
        rng=None,
    ):
        self.bond_length = bond_length
        self.radius = radius
        if bond_length < radius:
            raise ValueError(
                "Bond length should be greater than radius to prevent overlaps."
            )
        # Single RNG drives all walk randomness (angles, positions, bias,
        # volume-constraint sampling).
        if rng is None:
            rng = np.random.default_rng(seed + previous_count)
        self.rng = rng
        # Multiple ways to handle angles_sampler arg:
        if angles_sampler is None:
            self.angles = AnglesSampler(
                "uniform", {"low": np.pi / 2, "high": np.pi}, rng=self.rng
            )
        # Pass in a tupe or list of (low, high)
        elif isinstance(angles_sampler, (tuple, list)) and len(angles_sampler) == 2:
            self.angles = AnglesSampler(
                "uniform",
                {"low": angles_sampler[0], "high": angles_sampler[1]},
                rng=self.rng,
            )
        # Pass in a dict with supported kwargs
        elif isinstance(angles_sampler, dict):
            if "loc" in angles_sampler and "scale" in angles_sampler:
                self.angles = AnglesSampler("normal", angles_sampler, rng=self.rng)
            elif "low" in angles_sampler and "high" in angles_sampler:
                self.angles = AnglesSampler("uniform", angles_sampler, rng=self.rng)
            else:
                raise ValueError(
                    f"kwargs {angles_sampler} cannot be used to create an "
                    "AnglesSampler. Pass either {'loc': mean, 'scale': std} for "
                    "a normal distribution, or {'low': min, 'high': max} for a "
                    "uniform distribution."
                )
        # Pass in an array of choices
        elif isinstance(angles_sampler, np.ndarray):
            if angles_sampler.ndim == 1:
                kwargs = {"a": angles_sampler}
            elif angles_sampler.ndim == 2:
                kwargs = {"a": angles_sampler[0], "p": angles_sampler[1]}
            else:
                raise ValueError(
                    "Sampling angles from an array of choices is only supported for 1D and 2D arrays."
                )
            self.angles = AnglesSampler("choice", kwargs, rng=self.rng)
        # Pass in an AnglesSampler instance.
        elif isinstance(angles_sampler, AnglesSampler):
            self.angles = angles_sampler
            self.angles.rng = self.rng
        else:
            raise ValueError(
                f"{angles_sampler} is not a supported form to sample angles. "
                "See mbuild.path.points.AnglesSampler."
            )
        self.bead_name = bead_name
        self.initial_point, self.starting_from_site = _normalize_initial_point(
            initial_point
        )
        self.previous_count = previous_count
        self.include_compound = include_compound
        self.connectivity = connectivity
        self.seed = seed
        self.volume_constraint = volume_constraint
        self.termination = termination
        self.tolerance = tolerance
        self.bias = bias
        self.trial_batch_size = trial_batch_size
        self.chunk_size = chunk_size

        # State tracking
        self.count = 0
        self.init_count = 0
        self.attempts = 0
        self.start_time = None
        # PBC info for overlap checks; populated in hard_sphere_random_walk.
        # Defaults reproduce non-periodic behavior.
        self.pbc = np.array([False, False, False], dtype=np.bool_)
        self.box_lengths = np.array([np.inf, np.inf, np.inf], dtype=np.float32)

    def check_termination(self, path, coordinates, beads):
        """Examine and process termination if we have reached.

        This method also clears the path, state, and coordinate
        data structures from RW objects like Termination and Bias.

        Parameters
        ----------
        path : mbuild.path.Path, required
            The path used in the hard sphere random walk
        coordinates : np.ndarray (N, 3), required
            The up-to-date coordinates array from the last random walk chunk
        beads : np.ndarray (N, 3), required
            The up-to-date bead name array from the last random walk chunk
        """
        if self.termination.is_met(
            coordinates=coordinates[: self.count], names=beads[: self.count]
        ):
            if self.termination.success:
                logger.info("Random walk successful.")
                path.coordinates = coordinates[: self.count]
                path.beads = beads[: self.count]
            else:
                logger.warning("Random walk not successful.")
                logger.warning(self.termination.summarize())
                return True
            # RW is terminated and successful, update bond graph
            self.termination._clean()
            if self.bias:
                self.bias._clean()
            path._extend_bond_graph()
            if self.starting_from_site:
                # build from the given site instead of the last point
                path._connect_edges(
                    self.connectivity,
                    np.arange(self.previous_count, self.count),
                    self.initial_point,
                )
            else:  # build bond graph, and connect to last index in previous path coordinates
                path._connect_edges(
                    self.connectivity,
                    np.arange(self.previous_count, self.count),
                    self.previous_count,
                )
            # path._extend_beads(self.bead_name)
            return True
        return False


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
