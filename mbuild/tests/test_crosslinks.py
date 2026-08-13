import copy

import numpy as np
import pytest

from mbuild.exceptions import PathConvergenceError
from mbuild.path.build import Path
from mbuild.path.constraints import CuboidConstraint
from mbuild.path.crosslink import CrosslinkerGeometry, crosslink
from mbuild.tests.base_test import BaseTest


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
    def parallel_chains(self):
        # Two linear paths offset by two
        coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=np.float32)
        linear_path = Path(coords, bead_name="_A")
        linear_path._connect_edges("linear")
        linear_path.append_coordinates(coords + np.array([0, 2, 0]), bead_names="_A")
        linear_path.form_linear_bond_graph(np.arange(3, 6))
        return linear_path

    @pytest.fixture
    def periodic_path(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=np.float32)
        linear_path = Path(coords, bead_name="_A")
        linear_path._connect_edges("linear")
        linear_path.append_coordinates(coords + np.array([0, 8, 0]), bead_names="_A")
        linear_path.form_linear_bond_graph(np.arange(3, 6))
        return linear_path

    @pytest.fixture
    def periodic_cuboid(self):
        """10×10×10 fully periodic box centred at (1, 4, 0)."""
        return CuboidConstraint(
            10.0,
            10.0,
            10.0,
            center=[1.0, 4.0, 0.0],
            pbc=(True, True, True),
        )

    @pytest.fixture
    def one_site_crosslink(self):
        return CrosslinkerGeometry.single_site()

    @pytest.fixture
    def two_site_crosslink(self):
        return CrosslinkerGeometry.linear(connection_sites=[0, 1])  # default 2-site

    def test_single_site(self, parallel_chains, one_site_crosslink):
        path = parallel_chains
        clinker = one_site_crosslink
        crosslink(
            path, clinker, crosslink_bond_length=1, seed=3, backbone_bead_name="_A"
        )
        crosslink(
            path, clinker, crosslink_bond_length=1, seed=3, backbone_bead_name="_A"
        )
        crosslink(
            path, clinker, crosslink_bond_length=1, seed=3, backbone_bead_name="_A"
        )
        # Original 4 sites only
        assert len(path.coordinates) == 9
        assert len(path.bond_graph.edges()) == 10
        for e1, e2 in path.bond_graph.edges():
            distance = np.linalg.norm(path.coordinates[e1] - path.coordinates[e2])
            assert np.isclose(distance, 1)

    def test_single_site_extended(self, parallel_chains, one_site_crosslink):
        path = parallel_chains
        clinker = one_site_crosslink
        min_sep = 1
        crosslink(
            path,
            clinker,
            crosslink_bond_length=1.1,
            seed=3,
            backbone_bead_name="_A",
            max_backbone_degree=2,
            minimum_separation=min_sep,
            tolerance=0,
        )
        crosslink(
            path,
            clinker,
            crosslink_bond_length=1.1,
            seed=3,
            backbone_bead_name="_A",
            max_backbone_degree=2,
            minimum_separation=min_sep,
            tolerance=0,
        )
        crosslink(
            path,
            clinker,
            crosslink_bond_length=1.1,
            seed=3,
            backbone_bead_name="_A",
            max_backbone_degree=2,
            minimum_separation=min_sep,
            tolerance=0,
        )
        assert len(path.coordinates) == 9
        assert len(path.bond_graph.edges()) == 10
        print(path.bond_graph.edges())
        for e1, e2 in path.bond_graph.edges():
            distance = np.linalg.norm(path.coordinates[e1] - path.coordinates[e2])
            assert np.isclose(distance, 1) or np.isclose(distance, 1.1)

    def test_direct_crosslink_aperiodic(self, parallel_chains):
        """clinker=None: direct bond placed between chains 2 units apart,
        no periodic boundary conditions.  Validates that:
        - exactly one new bond is added (no new beads),
        - the bond length matches the target within tolerance.
        """
        path = parallel_chains
        n_beads_before = len(path.coordinates)
        n_edges_before = len(path.bond_graph.edges())
        bond_length = 2.1

        crosslink(
            path,
            crosslinker=None,
            crosslink_bond_length=bond_length,
            tolerance=0.1,
            seed=0,
            backbone_bead_name="_A",
        )

        assert len(path.coordinates) == n_beads_before
        assert len(path.bond_graph.edges()) == n_edges_before + 1

        # Check bond length tolerance
        for e1, e2 in path.bond_graph.edges():
            if abs(e1 - e2) < 3:
                continue
            distance = np.linalg.norm(path.coordinates[e1] - path.coordinates[e2])
            assert bond_length * 0.9 <= distance <= bond_length * 1.1, (
                f"Direct bond length {distance:.4f} outside tolerance of {bond_length}"
            )

    def test_direct_crosslink_periodic(self, periodic_path, periodic_cuboid):
        """clinker=None: direct bond placed across a periodic boundary.

        The two chains are 8 units apart in y inside a box of length 10,
        so the minimum-image distance is 2.  Without PBC the pair would
        be at distance 8 and the crosslink should *not* form.
        """
        path = periodic_path
        bond_length = 2.0
        n_beads_before = len(path.coordinates)
        n_edges_before = len(path.bond_graph.edges())

        # Fails with separation of 8, not 2 without pbc
        with pytest.raises(PathConvergenceError):
            crosslink(
                copy.deepcopy(path),
                crosslinker=None,
                crosslink_bond_length=bond_length,
                tolerance=0.1,
                seed=0,
                backbone_bead_name="_A",
                volume_constraint=None,
            )

        # With the periodic box the minimum-image distance is 2.
        crosslink(
            path,
            crosslinker=None,
            crosslink_bond_length=bond_length,
            tolerance=0.1,
            seed=0,
            backbone_bead_name="_A",
            volume_constraint=periodic_cuboid,
        )

        assert len(path.coordinates) == n_beads_before
        assert len(path.bond_graph.edges()) == n_edges_before + 1
        for e1, e2 in path.bond_graph.edges():
            test_distance = np.linalg.norm(path.coordinates[e1] - path.coordinates[e2])
            assert test_distance in (1.0, 8.0)

    def test_single_site_minimum_separation_respected(
        self, parallel_chains, one_site_crosslink
    ):
        """Placed crosslinker bead must be at least minimum_separation
        away from every pre-existing bead.
        """
        path = parallel_chains
        min_sep = 0.8

        crosslink(
            path,
            crosslinker=one_site_crosslink,
            crosslink_bond_length=1.0,
            seed=7,
            backbone_bead_name="_A",
            minimum_separation=min_sep,
            tolerance=0.1,
        )

        # Index of last added bead.
        new_bead_idx = len(path.coordinates) - 1
        new_bead_pos = path.coordinates[new_bead_idx]
        original_coords = path.coordinates[:new_bead_idx]

        distances = np.linalg.norm(original_coords - new_bead_pos, axis=1)
        assert np.all(distances >= min_sep), (
            f"Crosslinker bead placed too close to an existing bead: "
            f"min distance = {distances.min():.4f}, required >= {min_sep}"
        )

    def test_single_site_tolerance_enforced(self, parallel_chains, one_site_crosslink):
        """Tolerance=0 with a bond_length that has no exact solution raises
        PathConvergenceError; a slightly relaxed tolerance succeeds.

        The chains are exactly 2 units apart in y.  A target bond_length
        of 1.3 (each arm = 1.3) gives a total reach of 2.6 > 2, so the
        geometry *is* solvable — but only if tolerance > 0 allows the
        solver enough room.  With tolerance=0 the sphere-intersection
        circle collapses and no valid position can be found.
        """

        # tolerance=0: won't find correct bond length
        fail_path = copy.deepcopy(parallel_chains)
        with pytest.raises(PathConvergenceError):
            crosslink(
                fail_path,
                crosslinker=one_site_crosslink,
                crosslink_bond_length=0.5,
                tolerance=0.0,
                seed=1,
                backbone_bead_name="_A",
            )

        # tolerance=0.2: reaches bond lengths of one
        bond_length = 0.9
        crosslink(
            parallel_chains,
            crosslinker=one_site_crosslink,
            crosslink_bond_length=bond_length,
            tolerance=0.2,
            seed=1,
            backbone_bead_name="_A",
        )

        # Verify placed bond lengths sit within the relaxed tolerance window.
        new_bead_idx = len(parallel_chains.coordinates) - 1
        new_bead_pos = parallel_chains.coordinates[new_bead_idx]
        for nb in parallel_chains.bond_graph.neighbors(new_bead_idx):
            d = np.linalg.norm(parallel_chains.coordinates[nb] - new_bead_pos)
            assert bond_length * 0.8 <= d <= bond_length * 1.2, (
                f"Bond length {d:.4f} outside 20% tolerance of {bond_length}"
            )

    def test_crosslinker_wraps_into_box(self):
        """Build a periodic box and place beads on opposite side of periodic boundary."""
        box_len = 4.0
        half = box_len / 2.0
        near_face = half - 0.05  # 0.05 nm inside the +x face
        sep = 1.0  # inter-chain separation in y

        coords = np.array(
            [
                [near_face, 0.0, 0.0], 
                [-near_face, 0.0, 0.0],  
                [near_face, sep, 0.0], 
                [-near_face, sep, 0.0],  
            ],
            dtype=np.float32,
        )
        path = Path(coords, bead_name="_A")
        path.bond_graph.add_edge(0, 2)
        path.bond_graph.add_edge(1, 3)

        constraint = CuboidConstraint(
            box_len,
            box_len,
            box_len,
            center=[0.0, 0.0, 0.0],
            pbc=(True, True, True),
        )

        bulky_crosslinker = CrosslinkerGeometry.linear(
            n_sites=9,
            bond_length=0.4,
            bead_name="_R",
            connection_sites=0,
        )

        # MIC distance between bead 0 and bead 1 across the periodic x-boundary
        # is 2 * 0.05 = 0.1 nm — use that as the bond length so the
        # sphere-intersection solver is forced to place beads near that face.
        bond_length = 0.1
        n_beads_before = len(path.coordinates)
        crosslink(
            path,
            crosslinker=bulky_crosslinker,
            crosslink_bond_length=bond_length,
            tolerance=0.15,
            seed=42,
            backbone_bead_name="_A",
            minimum_intrachain_depth=1,
            volume_constraint=constraint,
        )

        assert len(path) == 13  # check crosslinker sites were added
        coords_after = path.coordinates
        pbc = np.array([True, True, True])

        violations = []
        for ax, (periodic, L) in enumerate(zip(pbc, [box_len, box_len, box_len])):
            if not periodic:
                continue
            lo, hi = -L / 2.0, L / 2.0
            out_lo = np.where(coords_after[:, ax] < lo - 1e-5)[0]
            out_hi = np.where(coords_after[:, ax] > hi + 1e-5)[0]
            for idx in np.concatenate([out_lo, out_hi]):
                violations.append(
                    f"  bead {idx:4d}  axis={ax}  pos={coords_after[idx, ax]:.6f}"
                    f"  box=[{lo:.3f}, {hi:.3f}]"
                )

        assert not violations, (
            f"{len(violations)} bead(s) found outside the periodic box after crosslinking:\n"
            + "\n".join(violations)
        )

        # Check all bond lengths
        new_bead_indices = np.arange(n_beads_before, len(path.coordinates))
        for new_idx in new_bead_indices:
            for nb in path.bond_graph.neighbors(new_idx):
                p_new = path.coordinates[new_idx].astype(np.float64)
                p_nb = path.coordinates[nb].astype(np.float64)
                # minimum-image distance so PBC bonds are measured correctly
                diff = p_new - p_nb
                for ax in range(3):
                    diff[ax] -= box_len * np.round(diff[ax] / box_len)
                d = np.linalg.norm(diff)
                # bonds to backbone anchors must respect the crosslink_bond_length
                if nb < n_beads_before:
                    assert d <= bond_length * (1 + 0.15) + 1e-4, (
                        f"Crosslinker→backbone bond too long: {d:.4f} "
                        f"(max {bond_length * 1.15:.4f})"
                    )

    def test_multiple_crosslinks(self, parallel_chains):
        crosslink(
            parallel_chains,
            None,
            crosslink_bond_length=2,
            n_crosslinks=3,
            minimum_intrachain_depth=2,
        )
        assert len(parallel_chains.bond_graph.edges()) == 7
        edges = set(list(parallel_chains.bond_graph.edges()))
        for edge in [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]:
            assert edge in edges
