import numpy as np
import pytest

from mbuild.exceptions import PathConvergenceError
from mbuild.path.build import Path
from mbuild.path.crosslink import CrosslinkerGeometry, crosslink
from mbuild.path.replace_sites import replace_sites
from mbuild.tests.base_test import BaseTest


class TestCLinks(BaseTest):
    @pytest.fixture
    def linear_path(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [3, 0, 0], [5, 0, 0]])
        return Path(coords)

    def test_clink_equilateral(self, linear_path):
        cl = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27,
            bead_name="_R",
            connection_sites=[0, 1],  # only sites 0 and 1 bond to backbone
        )
        outPath = crosslink(linear_path, crosslinker=cl, crosslink_bond_length=0.5)
        assert len(outPath) == 7
        assert len(outPath.bond_graph.edges()) == 5
        # Either edge formation is acceptable
        assert outPath.bond_graph.has_edge(0, 4) or outPath.bond_graph.has_edge(0, 5)
        assert outPath.bond_graph.has_edge(1, 5) or outPath.bond_graph.has_edge(1, 4)


class TestCrosslinkerGeometry(BaseTest):
    """Tests for CrosslinkerGeometry construction and factory methods."""

    def test_single_site_default(self):
        cl = CrosslinkerGeometry.single_site()
        assert cl.n_sites == 1
        assert cl.n_connections == 2
        assert cl.connection_sites == [0, 0]
        assert len(cl.internal_bonds) == 0

    def test_single_site_custom_connections(self):
        cl = CrosslinkerGeometry.single_site(bead_name="_X", n_connections=3)
        assert cl.n_sites == 1
        assert cl.n_connections == 3
        assert cl.beads[0] == "_X"

    def test_equilateral_triangle_geometry(self):
        cl = CrosslinkerGeometry.equilateral_triangle(bond_length=0.27)
        assert cl.n_sites == 3
        assert cl.n_connections == 3
        assert len(cl.internal_bonds) == 3
        # Verify equilateral: all edge lengths should be ~0.27
        for i, j in cl.bond_graph.edges():
            dist = np.linalg.norm(cl.coordinates[i] - cl.coordinates[j])
            assert abs(dist - 0.27) < 1e-4

    def test_equilateral_triangle_partial_connections(self):
        cl = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1]
        )
        assert cl.n_connections == 2
        assert cl.connection_sites == [0, 1]

    def test_linear_geometry(self):
        cl = CrosslinkerGeometry.linear(n_sites=3, bond_length=0.5)
        assert cl.n_sites == 3
        assert cl.n_connections == 2  # default: first and last
        assert len(cl.internal_bonds) == 2
        # Check bond lengths
        for i, j in cl.bond_graph.edges():
            dist = np.linalg.norm(cl.coordinates[i] - cl.coordinates[j])
            assert abs(dist - 0.5) < 1e-4

    def test_square_geometry(self):
        cl = CrosslinkerGeometry.square(bond_length=0.27)
        assert cl.n_sites == 4
        assert cl.n_connections == 4
        assert len(cl.internal_bonds) == 4
        for i, j in cl.bond_graph.edges():
            dist = np.linalg.norm(cl.coordinates[i] - cl.coordinates[j])
            assert abs(dist - 0.27) < 1e-4

    def test_tetrahedral_geometry(self):
        cl = CrosslinkerGeometry.tetrahedral(bond_length=0.27)
        assert cl.n_sites == 4
        assert cl.n_connections == 4
        assert len(cl.internal_bonds) == 6  # complete graph K4
        for i, j in cl.bond_graph.edges():
            dist = np.linalg.norm(cl.coordinates[i] - cl.coordinates[j])
            assert abs(dist - 0.27) < 1e-4

    def test_trigonal_bipyramidal_geometry(self):
        cl = CrosslinkerGeometry.trigonal_bipyramidal(bond_length=0.27)
        assert cl.n_sites == 5
        assert cl.n_connections == 5
        assert len(cl.internal_bonds) == 9

    def test_pentagon_geometry(self):
        cl = CrosslinkerGeometry.pentagon(bond_length=0.27)
        assert cl.n_sites == 5
        assert cl.n_connections == 5
        assert len(cl.internal_bonds) == 5
        for i, j in cl.bond_graph.edges():
            dist = np.linalg.norm(cl.coordinates[i] - cl.coordinates[j])
            assert abs(dist - 0.27) < 1e-4

    def test_from_edges(self):
        coords = [[0, 0, 0], [1, 0, 0], [0.5, 0.866, 0]]
        edges = [(0, 1), (1, 2)]
        cl = CrosslinkerGeometry.from_edges(
            coords, edges, bead_name="_R", connection_sites=[0, 2]
        )
        assert cl.n_sites == 3
        assert cl.n_connections == 2
        assert len(cl.internal_bonds) == 2
        assert cl.connection_sites == [0, 2]

    def test_from_path(self):
        coords = np.array([[0, 0, 0], [0.27, 0, 0]], dtype=np.float32)
        p = Path(coords, bead_name="_R")
        p._connect_edges("linear")
        cl = CrosslinkerGeometry.from_path(p, connection_sites=[0, 1])
        assert cl.n_sites == 2
        assert cl.n_connections == 2
        assert len(cl.internal_bonds) == 1

    def test_centered_at_origin(self):
        cl = CrosslinkerGeometry.equilateral_triangle(bond_length=0.27)
        centroid = np.mean(cl.coordinates, axis=0)
        assert np.allclose(centroid, [0, 0, 0], atol=1e-6)

    def test_invalid_connection_site_raises(self):
        coords = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)
        with pytest.raises(
            PathConvergenceError, match="connection site 5 is too large"
        ):
            CrosslinkerGeometry(
                coordinates=coords, bead_name="_R", connection_sites=[0, 5]
            )

    def test_unique_connection_sites(self):
        cl = CrosslinkerGeometry.single_site(n_connections=3)
        assert cl.unique_connection_sites == [0]
        cl2 = CrosslinkerGeometry.equilateral_triangle(connection_sites=[0, 1, 2])
        assert cl2.unique_connection_sites == [0, 1, 2]

    def test_copy(self):
        cl = CrosslinkerGeometry.equilateral_triangle(bond_length=0.27)
        cl2 = cl.copy()
        assert np.allclose(cl.coordinates, cl2.coordinates)
        assert cl.connection_sites == cl2.connection_sites
        # Mutate copy, original unchanged
        cl2.coordinates[0] = [99, 99, 99]
        assert not np.allclose(cl.coordinates, cl2.coordinates)


class TestCrosslink(BaseTest):
    """Tests for the crosslink function."""

    @pytest.fixture
    def linear_path(self):
        coords = np.array(
            [[0, 0, 0], [1, 0, 0], [3, 0, 0], [5, 0, 0]], dtype=np.float32
        )
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")
        return p

    @pytest.fixture
    def long_linear_path(self):
        n = 20
        coords = np.column_stack(
            (np.linspace(0, 5, n), np.zeros(n), np.zeros(n))
        ).astype(np.float32)
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")
        return p

    def test_crosslink_single_site_default(self, linear_path):
        out = crosslink(linear_path, crosslink_bond_length=3.0)
        # Original 4 + 1 crosslink = 5
        assert len(out.coordinates) == 5
        # 3 backbone edges + 2 crosslink edges = 5
        assert len(out.bond_graph.edges()) == 5

    def test_crosslink_equilateral_2conn(self, linear_path):
        cl = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1]
        )
        out = crosslink(linear_path, crosslinker=cl, crosslink_bond_length=3.0)
        # 4 backbone + 3 triangle sites = 7
        assert len(out.coordinates) == 7
        # 3 backbone + 3 internal + 2 connections = 8
        assert len(out.bond_graph.edges()) == 8

    def test_crosslink_equilateral_3conn(self, long_linear_path):
        cl = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1, 2]
        )
        out = crosslink(
            long_linear_path,
            crosslinker=cl,
            crosslink_bond_length=2.0,
            excluded_bond_depth=1,
        )
        # 20 backbone + 3 triangle = 23
        assert len(out.coordinates) == 23
        # 19 backbone edges + 3 internal + 3 connections = 25
        assert len(out.bond_graph.edges()) == 25

    def test_crosslink_linear_dimer(self, linear_path):
        cl = CrosslinkerGeometry.linear(
            n_sites=2, bond_length=0.27, connection_sites=[0, 1]
        )
        out = crosslink(linear_path, crosslinker=cl, crosslink_bond_length=3.0)
        # 4 + 2 = 6
        assert len(out.coordinates) == 6
        # 3 backbone + 1 internal + 2 connections = 6
        assert len(out.bond_graph.edges()) == 6

    def test_crosslink_no_backbone_raises(self):
        coords = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)
        p = Path(coords, bead_name="_B")  # not "_A"
        p._connect_edges("linear")
        with pytest.raises(
            PathConvergenceError,
            match="Could not find backbone beads matching crosslinker geometry",
        ):
            crosslink(p, backbone_name="_A", crosslink_bond_length=1.0)

    def test_crosslink_radius_too_small_raises(self, linear_path):
        with pytest.raises(PathConvergenceError):
            crosslink(linear_path, crosslink_bond_length=0.001)

    def test_crosslink_preserves_backbone_edges(self, linear_path):
        out = crosslink(linear_path, crosslink_bond_length=3.0)
        # Backbone edges 0-1, 1-2, 2-3 should still exist
        assert out.bond_graph.has_edge(0, 1)
        assert out.bond_graph.has_edge(1, 2)
        assert out.bond_graph.has_edge(2, 3)

    def test_crosslink_bead_names_correct(self, linear_path):
        out = crosslink(
            linear_path,
            crosslink_bond_length=3.0,
        )
        # Last bead should be the crosslinker
        assert out.beads[4] == "_R"
        # Backbone names unchanged
        for i in range(4):
            assert out.beads[i] == "_A"

    def test_crosslink_index_consistency(self, linear_path):
        out = crosslink(linear_path, crosslink_bond_length=3.0)
        # Every node in bond_graph should be a valid coordinate index
        for node in out.bond_graph.nodes():
            assert 0 <= node < len(out.coordinates)

    def test_clink_repeated_replaces(self, long_linear_path):
        triangle = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1], bead_name="_U"
        )
        crosslink(long_linear_path, triangle, crosslink_bond_length=1)
        crosslink(long_linear_path, triangle, crosslink_bond_length=1)
        assert len(long_linear_path) == 26
        assert sum(long_linear_path.beads == "_U") == 6
        assert sum(long_linear_path.beads == "_A") == 20
        assert long_linear_path.bond_graph.number_of_edges() == 29

        triangle = CrosslinkerGeometry.square(
            bond_length=0.27, connection_sites=[0, 1, 2, 3], bead_name="_U"
        )
        replace_sites(long_linear_path, triangle, bead_name="_A")
        assert len(long_linear_path) == 86
        assert sum(long_linear_path.beads == "_U") == 86

    def test_error_on_excess_degree(self, long_linear_path):
        triangle = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1], bead_name="_U"
        )
        crosslink(long_linear_path, triangle, crosslink_bond_length=5)
        with pytest.raises(ValueError, match="degree"):
            replace_sites(long_linear_path, triangle, bead_name="_A")

    def test_non_centroid_candidate(self, long_linear_path):
        triangle = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1], bead_name="_U"
        )
        crosslink(long_linear_path, triangle, crosslink_bond_length=5)
        assert long_linear_path

    # def test_high_density_crosslinking_overlaps_linear(self):
    #     # Make a gridded bond graph
    #     X, Y, Z = np.meshgrid(np.arange(5), np.arange(5), np.arange(5), indexing="ij")
    #     coords = np.column_stack((X.ravel(), Y.ravel(), Z.ravel())).astype(np.float32)
    #     path = Path(coordinates=coords)
    #     G = path.bond_graph
    #     for i in range(len(coords)):
    #         for j in range(i + 1, len(coords)):
    #             distance = np.linalg.norm(coords[i] - coords[j])
    #             if distance < 1.01:
    #                 G.add_edge(i, j)

    #     bond_length = np.sqrt(3) / 2  # half the diagonal bond length
    #     clink_coords = np.array([[0, 0, 0], [0, 0, bond_length / 2]], dtype=np.float32)
    #     clinkG = nx.Graph()
    #     clinkG.add_nodes_from([0, 1])
    #     clinkG.add_edge(0, 1)
    #     clink = CrosslinkerGeometry(
    #         clink_coords, clinkG, bead_name="_CR", connection_sites=[0, 0]
    #     )
    #     min_sep = 0.99 * (bond_length / 2)
    #     for i in range(4 * 4 * 4):
    #         try:
    #             crosslink(
    #                 path,
    #                 clink,
    #                 excluded_bond_depth=2,
    #                 minimum_separation=min_sep,
    #                 crosslink_bond_length=bond_length,
    #                 tolerance=0.01,
    #                 max_backbone_degree=8,
    #             )
    #         except PathConvergenceError:
    #             raise ValueError(f"Failed after {i} attempts to crosslink")

    #     # Find the non-connection crosslinker sites (site index 1 in the geometry)
    #     # These are "_CR" beads that are NOT the connection site (site 0)
    #     clink_nodes = [n for n in path.bond_graph.nodes() if path.beads[n] == "_CR"]
    #     # The connection site (index 0) has degree 2+1=3 (two backbone + one internal)
    #     # The internal site (index 1) has degree 1 (just the internal bond)
    #     internal_nodes = [n for n in clink_nodes if path.bond_graph.degree[n] == 1]

    #     assert len(internal_nodes) > 0, "No internal crosslinker sites found"

    #     # For each internal node, compute min distance to all non-bonded nodes
    #     for node in internal_nodes:
    #         bonded = set(path.bond_graph.neighbors(node)) | {node}
    #         node_pos = path.coordinates[node]

    #         other_positions = np.array(
    #             [
    #                 path.coordinates[i]
    #                 for i in path.bond_graph.nodes()
    #                 if i not in bonded
    #             ],
    #             dtype=np.float32,
    #         )

    #         distances = np.linalg.norm(other_positions - node_pos, axis=1)
    #         min_dist = float(np.min(distances))

    #         assert min_dist >= min_sep, (
    #             f"Internal crosslinker node {node} has overlap: "
    #             f"min distance {min_dist:.4f} < minimum_separation {min_sep}"
    #         )

    # def test_high_density_crosslinking_overlaps_eqtriangle(self):
    #     # Make a gridded bond graph
    #     X, Y, Z = np.meshgrid(np.arange(5), np.arange(5), np.arange(5), indexing="ij")
    #     coords = np.column_stack((X.ravel(), Y.ravel(), Z.ravel())).astype(np.float32)
    #     path = Path(coordinates=coords)
    #     G = path.bond_graph
    #     for i in range(len(coords)):
    #         for j in range(i + 1, len(coords)):
    #             distance = np.linalg.norm(coords[i] - coords[j])
    #             if distance < 1.01:
    #                 G.add_edge(i, j)

    #     bond_length = np.sqrt(3) / 3  # 1/3 the diagonal bond length
    #     clink = CrosslinkerGeometry.equilateral_triangle(
    #         bond_length=bond_length, connection_sites=[0, 1]
    #     )
    #     min_sep = 0.1
    #     for i in range(4 * 4 * 4):
    #         crosslink(
    #             path,
    #             clink,
    #             excluded_bond_depth=2,
    #             minimum_separation=min_sep,
    #             crosslink_bond_length=bond_length,
    #             tolerance=1,
    #             max_backbone_degree=8,
    #         )

    #     # Find the non-connection crosslinker sites (site index 1 in the geometry)
    #     # These are "_CR" beads that are NOT the connection site (site 0)
    #     clink_nodes = [n for n in path.bond_graph.nodes() if path.beads[n] == "_R"]
    #     # The connection site (index 0) has degree 2+1=3 (two backbone + one internal)
    #     # The internal site (index 1) has degree 1 (just the internal bond)
    #     internal_nodes = [n for n in clink_nodes if path.bond_graph.degree[n] == 2]

    #     assert len(internal_nodes) > 0, "No internal crosslinker sites found"

    #     # For each internal node, compute min distance to all non-bonded nodes
    #     for node in internal_nodes:
    #         bonded = set(path.bond_graph.neighbors(node)) | {node}
    #         node_pos = path.coordinates[node]

    #         other_positions = np.array(
    #             [
    #                 path.coordinates[i]
    #                 for i in path.bond_graph.nodes()
    #                 if i not in bonded
    #             ],
    #             dtype=np.float32,
    #         )

    #         distances = np.linalg.norm(other_positions - node_pos, axis=1)
    #         min_dist = float(np.min(distances))

    #         assert min_dist >= min_sep, (
    #             f"Internal crosslinker node {node} has overlap: "
    #             f"min distance {min_dist:.4f} < minimum_separation {min_sep}"
    #         )

    # def test_failing_overlap_radius(self, long_linear_path):
    #     clink_coords = np.array([[0, 0, 0], [0, 0, 0.5]], dtype=np.float32)
    #     clinkG = nx.Graph()
    #     clinkG.add_nodes_from([0, 1])
    #     clinkG.add_edge(0, 1)
    #     linear_clink = CrosslinkerGeometry(
    #         clink_coords, clinkG, connection_sites=[0, 1]
    #     )
    #     with pytest.raises(PathConvergenceError):
    #         crosslink(
    #             long_linear_path,
    #             linear_clink,
    #             crosslink_bond_length=0.5,
    #             minimum_separation=0.2,
    #         )

    def test_multiple_connection_sites(self):
        coordinates = np.array([[0, 0, 0], [1, 1, 0], [1, 0, 0], [0, 1, 0]])
        path = Path(coordinates, bead_name="_B")
        path.bond_graph.add_edges_from([[0, 2], [1, 3]])
        cl = CrosslinkerGeometry(
            bead_name="_CROSS",
            connection_sites=[0, 0],
            coordinates=np.array([[0, 0, 0]]),
        )
        crosslink(
            path,
            crosslinker=cl,
            backbone_name=(("_B", "_B"), ("_B", "_B")),
            crosslink_bond_length=np.sqrt(2) / 2,
            tolerance=0.01,
            seed=1,
        )
        assert len(path) == 5  # 4 _B and one _CR
        assert len(path.bond_graph.edges) == 6  # 4 _B-_CROSS bonds
        assert np.allclose(path.coordinates[-1], np.array([0.5, 0.5, 0]))  # CROSS Bead
        assert (
            np.linalg.norm(path.coordinates[-1] - path.coordinates[-2])
            == np.sqrt(2) / 2
        )

    def test_crosslink_no_overlaps_sparse(self):
        """Sparse system: crosslinks should place without overlaps."""
        # 4 chains of 5 beads in a large box
        coords = []
        for i in range(4):
            for j in range(5):
                coords.append([i * 3.0, j * 0.5, 0.0])
        coords = np.array(coords, dtype=np.float32)
        path = Path(coordinates=coords, bead_name="A")
        # Connect each chain linearly
        for chain in range(4):
            for j in range(4):
                path.bond_graph.add_edge(chain * 5 + j, chain * 5 + j + 1)

        radius = 0.2
        cl = CrosslinkerGeometry.single_site(bead_name="R")
        crosslink(
            path,
            crosslinker=cl,
            n_crosslinks=2,
            crosslink_bond_length=1.6,
            minimum_separation=radius,
            tolerance=0.1,
            seed=42,
        )

        # Check no overlaps exist (excluding bonded pairs)
        n_total = len(path.coordinates)
        for i in range(n_total):
            for j in range(i + 1, n_total):
                if path.bond_graph.has_edge(i, j):
                    continue
                dist = np.linalg.norm(path.coordinates[i] - path.coordinates[j])
                assert dist >= radius * 0.99, (
                    f"Overlap between {i} and {j}: dist={dist:.4f} < {radius}"
                )

    def test_crosslink_dense_system_raises_or_limited(self):
        """Very dense system: should raise or place fewer crosslinks."""
        # Tightly packed beads where minimum_separation can't be satisfied
        n = 10
        coords = np.column_stack(
            [np.linspace(0, 2, n), np.zeros(n), np.zeros(n)]
        ).astype(np.float32)
        path = Path(coordinates=coords, bead_name="A")
        path._connect_edges("linear")

        cl = CrosslinkerGeometry.single_site(bead_name="R")
        # Bond length ~0.22, separation 0.2, beads are 0.22 apart
        # Very tight: should struggle or fail
        with pytest.raises(PathConvergenceError):
            crosslink(
                path,
                crosslinker=cl,
                n_crosslinks=5,  # way too many for 10 beads
                crosslink_bond_length=0.3,
                minimum_separation=0.2,
                tolerance=0.01,
                n_rotation_samples=20,
                seed=42,
            )

    def test_crosslink_respects_minimum_separation(self):
        """Placed crosslinker must be >= minimum_separation from non-bonded beads."""
        # Two parallel chains close enough to crosslink
        coords = []
        # Chain 0: y=0
        for i in range(10):
            coords.append([i * 0.4, 0.0, 0.0])
        # Chain 1: y=1.0
        for i in range(10):
            coords.append([i * 0.4, 1.0, 0.0])
        coords = np.array(coords, dtype=np.float32)
        path = Path(coordinates=coords, bead_name="A")
        # Connect each chain
        for j in range(9):
            path.bond_graph.add_edge(j, j + 1)
            path.bond_graph.add_edge(10 + j, 10 + j + 1)

        min_sep = 0.15
        cl = CrosslinkerGeometry.single_site(bead_name="R")
        crosslink(
            path,
            crosslinker=cl,
            n_crosslinks=2,
            crosslink_bond_length=0.6,
            minimum_separation=min_sep,
            tolerance=0.2,
            n_rotation_samples=72,
            seed=42,
        )

        # Verify minimum separation for all non-bonded pairs
        n_total = len(path.coordinates)
        for i in range(n_total):
            for j in range(i + 1, n_total):
                if path.bond_graph.has_edge(i, j):
                    continue
                dist = np.linalg.norm(path.coordinates[i] - path.coordinates[j])
                assert dist >= min_sep * 0.99, (
                    f"Overlap: nodes {i},{j} dist={dist:.4f} < min_sep={min_sep}"
                )


# =============================================================================
# Sub-function unit tests
# =============================================================================


class TestCrosslinkSubFunctions(BaseTest):
    """Unit tests for internal crosslink helper functions."""

    def test_find_equidistant_point_single_bead(self):
        """Single bead: point should be exactly target_distance away."""
        from mbuild.path.crosslink import _find_equidistant_point

        bead = np.array([[1.0, 2.0, 3.0]])
        target = 0.5
        point = _find_equidistant_point(bead, target)
        dist = np.linalg.norm(point - bead[0])
        assert np.isclose(dist, target, atol=1e-6)

    def test_find_equidistant_point_two_beads(self):
        """Two beads: point should be equidistant from both."""
        from mbuild.path.crosslink import _find_equidistant_point

        beads = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        target = 0.7
        point = _find_equidistant_point(beads, target)
        d0 = np.linalg.norm(point - beads[0])
        d1 = np.linalg.norm(point - beads[1])
        assert np.isclose(d0, target, atol=1e-4)
        assert np.isclose(d1, target, atol=1e-4)

    def test_find_equidistant_point_symmetric_square(self):
        """Four symmetric beads: centroid should be the answer when r matches."""
        from mbuild.path.crosslink import _find_equidistant_point

        # Square at distance sqrt(2)/2 from centroid
        beads = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        centroid = beads.mean(axis=0)
        r = np.linalg.norm(beads[0] - centroid)
        point = _find_equidistant_point(beads, r)
        assert np.allclose(point, centroid, atol=1e-4)

    def test_can_reach_all_beads_feasible(self):
        """Two beads within 2*r should be feasible."""
        from mbuild.path.crosslink import _can_reach_all_beads

        beads = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
        pbc = np.array([False, False, False])
        box = np.array([0.0, 0.0, 0.0])
        feasible, target = _can_reach_all_beads(beads, 0.5, 0.1, pbc, box)
        assert feasible
        assert target is not None
        # Verify distances
        d0 = np.linalg.norm(target - beads[0])
        d1 = np.linalg.norm(target - beads[1])
        assert np.isclose(d0, 0.5, atol=0.05)
        assert np.isclose(d1, 0.5, atol=0.05)

    def test_can_reach_all_beads_infeasible(self):
        """Two beads too far apart should be infeasible."""
        from mbuild.path.crosslink import _can_reach_all_beads

        beads = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
        pbc = np.array([False, False, False])
        box = np.array([0.0, 0.0, 0.0])
        feasible, target = _can_reach_all_beads(beads, 0.5, 0.1, pbc, box)
        assert not feasible
        assert target is None

    def test_has_overlaps_true(self):
        """Placed coords too close to existing should flag overlap."""
        from mbuild.path.crosslink import _has_overlaps

        placed = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        existing = np.array([[0.05, 0.0, 0.0]], dtype=np.float64)
        pbc = np.array([False, False, False])
        box = np.array([0.0, 0.0, 0.0])
        assert _has_overlaps(placed, existing, [], 0.1, pbc, box)

    def test_has_overlaps_false(self):
        """Placed coords far from existing should not flag overlap."""
        from mbuild.path.crosslink import _has_overlaps

        placed = np.array([[5.0, 5.0, 5.0]], dtype=np.float64)
        existing = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        pbc = np.array([False, False, False])
        box = np.array([0.0, 0.0, 0.0])
        assert not _has_overlaps(placed, existing, [], 0.1, pbc, box)

    def test_has_overlaps_excluded_indices(self):
        """Excluded indices should not trigger overlap."""
        from mbuild.path.crosslink import _has_overlaps

        placed = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        existing = np.array([[0.05, 0.0, 0.0]], dtype=np.float64)
        pbc = np.array([False, False, False])
        box = np.array([0.0, 0.0, 0.0])
        # Exclude index 0 (the only existing coord)
        assert not _has_overlaps(placed, existing, [0], 0.1, pbc, box)

    def test_parse_backbone_spec_string(self):
        """String spec should replicate to all sites."""
        from mbuild.path.crosslink import _parse_backbone_spec

        result = _parse_backbone_spec("A", [0, 0, 1])
        assert result == ["A", "A", "A"]

    def test_parse_backbone_spec_list(self):
        """List spec must match connection site count."""
        from mbuild.path.crosslink import _parse_backbone_spec

        result = _parse_backbone_spec(["A", "B"], [0, 1])
        assert result == ["A", "B"]

    def test_parse_backbone_spec_mismatch_raises(self):
        """Mismatched length should raise."""
        from mbuild.path.crosslink import _parse_backbone_spec

        with pytest.raises(PathConvergenceError):
            _parse_backbone_spec(["A", "B", "C"], [0, 1])

    def test_bead_matches_spec_with_underscore(self):
        """Bead matching should handle underscore prefix."""
        from mbuild.path.crosslink import _bead_matches_spec

        assert _bead_matches_spec("_A", "A")
        assert _bead_matches_spec("_A", "_A")
        assert _bead_matches_spec("A", "_A")
        assert not _bead_matches_spec("_B", "A")

    def test_bead_matches_spec_tuple(self):
        """Tuple spec should match if bead matches any element."""
        from mbuild.path.crosslink import _bead_matches_spec

        assert _bead_matches_spec("_A", ("A", "B"))
        assert _bead_matches_spec("_B", ("A", "B"))
        assert not _bead_matches_spec("_C", ("A", "B"))

    def test_rotation_about_axis_identity(self):
        """Zero angle rotation should be identity."""
        from mbuild.path.crosslink import _rotation_about_axis

        R = _rotation_about_axis(np.array([0, 0, 1.0]), 0.0)
        assert np.allclose(R, np.eye(3), atol=1e-10)

    def test_rotation_about_axis_90deg(self):
        """90° about z should rotate x to y."""
        from mbuild.path.crosslink import _rotation_about_axis

        R = _rotation_about_axis(np.array([0, 0, 1.0]), np.pi / 2)
        x_hat = np.array([1, 0, 0], dtype=np.float64)
        rotated = R @ x_hat
        assert np.allclose(rotated, [0, 1, 0], atol=1e-10)

    def test_rotation_matrix_between_vectors(self):
        """Should rotate a to b."""
        from mbuild.path.crosslink import _rotation_matrix_between_vectors

        a = np.array([1, 0, 0], dtype=np.float64)
        b = np.array([0, 1, 0], dtype=np.float64)
        R = _rotation_matrix_between_vectors(a, b)
        result = R @ a
        assert np.allclose(result, b, atol=1e-10)

    def test_rotation_matrix_between_same_vectors(self):
        """Same vector should give identity."""
        from mbuild.path.crosslink import _rotation_matrix_between_vectors

        a = np.array([0, 0, 1], dtype=np.float64)
        R = _rotation_matrix_between_vectors(a, a)
        assert np.allclose(R, np.eye(3), atol=1e-10)

    def test_pbc_delta_no_pbc(self):
        """Without PBC, delta is just a - b."""
        from mbuild.path.crosslink import _pbc_delta

        a = np.array([1.0, 2.0, 3.0])
        b = np.array([0.5, 1.0, 1.5])
        pbc = np.array([False, False, False])
        box = np.array([10.0, 10.0, 10.0])
        delta = _pbc_delta(a, b, box, pbc)
        assert np.allclose(delta, a - b)

    def test_pbc_delta_with_wrapping(self):
        """With PBC, delta should use minimum image."""
        from mbuild.path.crosslink import _pbc_delta

        a = np.array([9.5, 0.0, 0.0])
        b = np.array([0.5, 0.0, 0.0])
        pbc = np.array([True, False, False])
        box = np.array([10.0, 10.0, 10.0])
        delta = _pbc_delta(a, b, box, pbc)
        # Minimum image: 9.5 - 0.5 = 9.0, wrapped to -1.0
        assert np.isclose(delta[0], -1.0, atol=1e-10)


class TestReplaceSites(BaseTest):
    """Tests for the replace_sites function."""

    @pytest.fixture
    def chain_with_crosslinks(self):
        """A linear backbone with two single-site crosslinks."""
        n = 10
        coords = np.column_stack(
            (np.linspace(0, 5, n), np.zeros(n), np.zeros(n))
        ).astype(np.float32)
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")
        # Add two crosslinks
        p = crosslink(p, crosslink_bond_length=3.0, seed=42)
        p = crosslink(p, crosslink_bond_length=3.0, seed=99)
        return p

    def test_replace_single_site_with_triangle(self):
        # Build a simple path with a degree-2 node in the middle
        coords = np.array(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]],
            dtype=np.float32,
        )
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")
        # Change middle node to _R
        p.beads[2] = "_R"

        triangle = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1]
        )
        out = replace_sites(p, triangle, bead_name="_R")
        # 5 - 1 removed + 3 added = 7
        assert len(out.coordinates) == 7
        # Check all node indices are valid
        for node in out.bond_graph.nodes():
            assert 0 <= node < len(out.coordinates)

    def test_replace_by_index(self):
        coords = np.array(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=np.float32
        )
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")

        linear_repl = CrosslinkerGeometry.linear(
            n_sites=2, bond_length=0.3, connection_sites=[0, 1]
        )
        # Replace node 1 (degree 2)
        out = replace_sites(p, linear_repl, sites=[1])
        # 4 - 1 + 2 = 5
        assert len(out.coordinates) == 5
        assert (
            len(out.bond_graph.edges()) == 4
        )  # 2 kept + 1 internal + 1 external? let's just check > 0
        assert len(out.bond_graph.edges()) >= 4

    def test_replace_preserves_nonreplaced_edges(self):
        coords = np.array(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]],
            dtype=np.float32,
        )
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")
        p.beads[2] = "_R"

        triangle = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1]
        )
        out = replace_sites(p, triangle, bead_name="_R")

        # Original edge 0-1 and 3-4 should still exist in some form
        # (nodes may be renumbered but connectivity is preserved)
        # Check total number of edges: 2 kept backbone + 3 internal + 2 external = 7
        # Actually: edges 0-1, 1-2, 2-3, 3-4. Remove node 2. Kept: 0-1, 3-4.
        # Triangle: 3 internal + 2 connections to (what was 1 and 3)
        # Total: 2 + 3 + 2 = 7
        assert len(out.bond_graph.edges()) == 7

    def test_replace_multiple_sites(self):
        coords = np.array(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0], [5, 0, 0]],
            dtype=np.float32,
        )
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")
        p.beads[1] = "_R"
        p.beads[4] = "_R"

        linear_repl = CrosslinkerGeometry.linear(
            n_sites=2, bond_length=0.3, connection_sites=[0, 1]
        )
        out = replace_sites(p, linear_repl, bead_name="_R")
        # 6 - 2 + 4 = 8
        assert len(out.coordinates) == 8

    def test_replace_no_sites_noop(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=np.float32)
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")

        triangle = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1]
        )
        # No "_R" beads exist
        out = replace_sites(p, triangle, bead_name="_R")
        assert len(out.coordinates) == 3
        assert len(out.bond_graph.edges()) == 2

    def test_replace_neither_sites_nor_bead_raises(self):
        coords = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)
        p = Path(coords, bead_name="_A")
        triangle = CrosslinkerGeometry.equilateral_triangle(bond_length=0.27)
        with pytest.raises(ValueError, match="Must specify"):
            replace_sites(p, triangle)

    def test_replace_both_sites_and_bead_raises(self):
        coords = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)
        p = Path(coords, bead_name="_A")
        triangle = CrosslinkerGeometry.equilateral_triangle(bond_length=0.27)
        with pytest.raises(ValueError, match="Cannot specify both"):
            replace_sites(p, triangle, sites=[0], bead_name="_A")

    def test_replace_renumbers_bond_graph(self):
        coords = np.array(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=np.float32
        )
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")
        p.beads[2] = "_R"

        linear_repl = CrosslinkerGeometry.linear(
            n_sites=2, bond_length=0.3, connection_sites=[0, 1]
        )
        out = replace_sites(p, linear_repl, bead_name="_R")

        # All node indices should be contiguous 0..N-1
        nodes = sorted(out.bond_graph.nodes())
        assert nodes == list(range(len(out.coordinates)))

    def test_replace_via_path_method(self):
        coords = np.array(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=np.float32
        )
        p = Path(coords, bead_name="_A")
        p._connect_edges("linear")
        p.beads[1] = "_R"

        linear_repl = CrosslinkerGeometry.linear(
            n_sites=2, bond_length=0.3, connection_sites=[0, 1]
        )
        out = replace_sites(p, linear_repl, bead_name="_R")
        assert len(out.coordinates) == 5

    def test_crosslink_then_replace(self, chain_with_crosslinks):
        """Replace single-site crosslinks with triangles."""
        p = chain_with_crosslinks
        n_R = sum(1 for b in p.beads if b == "_R")
        # Each _R replaced: remove 1, add 3 -> net +2 per replacement
        expected_len = len(p.coordinates) - n_R + n_R * 3
        assert n_R >= 1

        triangle = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.27, connection_sites=[0, 1]
        )
        out = replace_sites(p, triangle, bead_name="_R")
        assert len(out.coordinates) == expected_len

        # No _R beads remain (they've been replaced with triangle beads)
        # Triangle beads are also "_R" so they should still exist
        # But original single-site _R nodes are gone
        for node in out.bond_graph.nodes():
            assert 0 <= node < len(out.coordinates)

    def test_replace_site(self):
        coords = np.column_stack(
            (np.arange(0, 5), np.zeros(5, dtype=int), np.zeros(5, dtype=int))
        )
        path = Path(coords)
        path._connect_edges("linear")
        assert len(path.bond_graph.edges()) == 4

        triangle = CrosslinkerGeometry.equilateral_triangle(
            bond_length=0.4, connection_sites=[0, 1]
        )
        path = replace_sites(path, triangle, 1)
        assert len(path) == 7
