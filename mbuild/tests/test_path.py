import numpy as np
import pytest

import mbuild as mb
from mbuild.exceptions import PathConvergenceError
from mbuild.path.build import (
    Path,
    cyclic,
    hard_sphere_random_walk,
    helix,
    knot,
    lamellar,
    spiral_2D,
    straight_line,
    zigzag,
)
from mbuild.path.constraints import (
    CuboidConstraint,
    CylinderConstraint,
    SphereConstraint,
)
from mbuild.path.namers import CyclicNamer, RandomNamer
from mbuild.path.path_utils import (
    check_path,
    local_density,
    target_density,
    target_sq_distances,
)
from mbuild.path.points import AnglesSampler
from mbuild.path.termination import (
    NumAttempts,
    NumSites,
    Termination,
)
from mbuild.tests.base_test import BaseTest
from mbuild.utils.geometry import bounding_box


class TestPaths(BaseTest):
    def test_from_coordinates(self):
        coords = np.random.uniform(-5, 5, size=(20, 3))
        path = Path(coordinates=coords, bead_name="X")
        assert np.array_equal(coords, path.coordinates)
        assert path.bond_graph.number_of_nodes() == len(coords)
        assert len(path.beads) == len(coords)

    def test_from_compound(self):
        compound = mb.Compound()
        last_site = None
        for i in range(10):
            if i % 2 == 0:
                this_site = mb.Compound(name="A", pos=(i, 0, 0))
            else:
                this_site = mb.Compound(name="B", pos=(i, 0, 0))
            compound.add(this_site)
            if last_site:
                compound.add_bond([this_site, last_site])
            last_site = this_site

        path = Path.from_compound(compound)
        assert np.array_equal(path.coordinates, compound.xyz)
        for idx, node in enumerate(path.bond_graph.nodes(data=True)):
            assert node[1]["name"] == compound[idx].name
        for edge1, edge2 in zip(path.bond_graph.edges(), compound.bond_graph.edges()):
            assert edge1[0] == compound.get_child_indices(edge2[0])[0]

    def test_straight_line(self):
        path = Path()  # works with empty path
        straight_line(path=path, spacing=0.20, N=5, direction=(1, 0, 0))
        assert len(path.coordinates) == 5
        assert path.bond_graph.number_of_edges() == 4
        # 5 sites = 4 bonds at 0.20 each
        assert np.allclose(
            np.linalg.norm(path.coordinates[-1] - path.coordinates[0]), 0.80
        )
        for edge in path.bond_graph.edges(data=True):
            assert np.allclose(edge[2]["direction"], np.array([1, 0, 0]))

    def test_initial_point(self):
        """Every builder places its first site at initial_point."""
        initial_point = (1.5, -2.0, 3.0)
        builders = [
            (straight_line, dict(spacing=0.25, N=10)),
            (cyclic, dict(spacing=0.25, N=20)),
            (knot, dict(spacing=0.25, N=40, m=3)),
            (spiral_2D, dict(N=30, a=0.5, b=0.1, spacing=0.25)),
            (zigzag, dict(N=20, spacing=0.25)),
            (helix, dict(N=30, radius=0.5, rise=0.1, twist=30)),
        ]
        for builder, kwargs in builders:
            at_origin = builder(**kwargs).coordinates
            shifted = builder(initial_point=initial_point, **kwargs).coordinates
            assert np.allclose(shifted[0], initial_point)
            # Only translated; the shape of the path is unchanged
            assert np.allclose(shifted - shifted[0], at_origin - at_origin[0])

    def test_initial_point_appends_to_path(self):
        """Builders place their segment at initial_point within an existing path."""
        path = straight_line(spacing=0.25, N=5)
        straight_line(path=path, spacing=0.25, N=5, initial_point=(4.0, 0.0, 0.0))
        assert len(path.coordinates) == 10
        assert np.allclose(path.coordinates[0], (0, 0, 0))
        assert np.allclose(path.coordinates[5], (4.0, 0.0, 0.0))
        # The two segments are bonded separately, not to each other
        assert path.bond_graph.number_of_edges() == 8

    def test_knot_not_closed(self):
        path = knot(spacing=0.25, N=50, m=3, closed=False)
        assert len(path.coordinates) == 50
        assert path.bond_graph.number_of_edges() == 49

    def test_cyclic_parameters(self):
        path = Path()
        cyclic(path=path, spacing=1, N=20)
        # C = 2*pi*r
        radius = 20 * 1 / (2 * np.pi)
        # Check that computed radius matches expected
        actual_radius = np.linalg.norm(path.coordinates[0])
        assert np.allclose(actual_radius, radius, atol=1e-2)

        path2 = Path()
        cyclic(path=path2, spacing=1, radius=10 / np.pi, N=None)
        assert len(path2.coordinates) == 20

        path3 = Path()
        cyclic(path=path3, N=20, radius=10 / np.pi, spacing=None)
        # Check spacing
        dist = np.linalg.norm(path3.coordinates[1] - path3.coordinates[0])
        assert np.allclose(dist, 1.0, atol=1e-2)

    def test_cyclic_bonding(self):
        path = Path()
        cyclic(path=path, spacing=1, N=20)
        assert path.bond_graph.number_of_edges() == 20
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles

    def test_knot(self):
        path = Path()
        knot(path=path, spacing=0.25, N=50, m=3)
        assert path.bond_graph.number_of_edges() == 50
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles

    def test_knot_bad_arg(self):
        path = Path()
        with pytest.raises(ValueError):
            knot(path=path, spacing=0.25, N=50, m=2)

    def test_spiral(self):
        path = Path()
        spiral_2D(path=path, N=50, a=0.5, b=2, spacing=0.25)
        assert path.bond_graph.number_of_edges() == 49
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles - 1

    def test_lamellar(self):
        path = Path()
        lamellar(
            path=path,
            spacing=0.25,
            num_layers=3,
            layer_separation=1.0,
            layer_length=3.0,
            num_stacks=3,
            stack_separation=1.0,
        )
        assert path.bond_graph.number_of_edges() == len(path.coordinates) - 1
        compound = path.to_compound()
        Lx, Ly, Lz = compound.get_boundingbox().lengths
        # The params used here should create a cubic-like lamellar structure
        # Y-direction will be slightly larger because of curves between layers
        assert np.allclose(Lx, Lz, atol=0.1)  # stacking and layering directions
        assert Ly > Lx

    def test_lamellar_direction(self):
        path_left_to_right = Path()
        lamellar(
            path=path_left_to_right,
            spacing=0.25,
            num_layers=3,
            layer_separation=1.0,
            layer_length=3.0,
            num_stacks=3,
            stack_separation=1.0,
            initial_point=(0, 0, 0),
        )

        path_right_to_left = Path()
        lamellar(
            path=path_right_to_left,
            spacing=0.25,
            num_layers=3,
            layer_separation=1.0,
            layer_length=3.0,
            num_stacks=3,
            stack_separation=1.0,
            initial_point=(0, 0, 0),
            left_to_right=False,
        )

        assert np.array_equal(
            path_left_to_right.coordinates[0], path_right_to_left.coordinates[0]
        )
        assert path_right_to_left.coordinates[1][1] < 0
        assert path_left_to_right.coordinates[1][1] > 0

    def test_lamellar_initial_point(self):
        path = Path()
        lamellar(
            path=path,
            spacing=0.25,
            num_layers=3,
            layer_separation=1.0,
            layer_length=3.0,
            num_stacks=3,
            stack_separation=1.0,
            initial_point=(1, 1, 1),
        )
        assert np.array_equal(path.coordinates[0], np.array([1, 1, 1]))
        assert np.allclose(path.coordinates[-1][2], 3.0, atol=0.5)

    def test_helix(self):
        path = Path()
        helix(path=path, N=50, radius=2.0, rise=0.5, twist=30)
        assert len(path.coordinates) == 50
        assert path.bond_graph.number_of_edges() == 49
        # Check that all points are roughly at radius distance from z-axis
        radii = np.sqrt(path.coordinates[:, 0] ** 2 + path.coordinates[:, 1] ** 2)
        assert np.allclose(radii, 2.0, atol=1e-5)

    def test_zigzag(self):
        path = Path()
        zigzag(path=path, N=20, spacing=1.0, angle_deg=120.0, sites_per_segment=5)
        assert len(path.coordinates) == 20
        assert path.bond_graph.number_of_edges() == 19

    def test_straight_line_cyclic_namer_alternating(self):
        path = Path()
        straight_line(path=path, spacing=0.2, N=6, bead_name=CyclicNamer(["_A", "_B"]))
        assert list(path.beads) == ["_A", "_B", "_A", "_B", "_A", "_B"]

    def test_straight_line_cyclic_namer_blocks(self):
        path = Path()
        straight_line(
            path=path, spacing=0.2, N=6, bead_name=CyclicNamer([("_A", 3), ("_B", 3)])
        )
        assert list(path.beads) == ["_A", "_A", "_A", "_B", "_B", "_B"]

    def test_cyclic_path_cyclic_namer(self):
        path = Path()
        cyclic(path=path, spacing=1, N=6, bead_name=CyclicNamer(["_A", "_B"]))
        assert list(path.beads) == ["_A", "_B", "_A", "_B", "_A", "_B"]

    def test_helix_cyclic_namer(self):
        path = Path()
        helix(
            path=path,
            N=4,
            radius=1.0,
            rise=0.5,
            twist=90,
            bead_name=CyclicNamer(["_A", "_B"]),
        )
        assert list(path.beads) == ["_A", "_B", "_A", "_B"]

    def test_zigzag_cyclic_namer(self):
        path = Path()
        zigzag(
            path=path,
            N=4,
            spacing=1.0,
            sites_per_segment=2,
            bead_name=CyclicNamer(["_A", "_B"]),
        )
        assert list(path.beads) == ["_A", "_B", "_A", "_B"]

    def test_namer_beads_in_bond_graph(self):
        path = Path()
        straight_line(path=path, spacing=0.2, N=4, bead_name=CyclicNamer(["_A", "_B"]))
        node_names = [d["name"] for _, d in path.bond_graph.nodes(data=True)]
        assert node_names == ["_A", "_B", "_A", "_B"]

    def test_string_bead_name_still_works(self):
        path = Path()
        straight_line(path=path, spacing=0.2, N=4, bead_name="_X")
        assert list(path.beads) == ["_X", "_X", "_X", "_X"]

    def test_add_paths(self):
        coords = np.random.uniform(-5, 5, size=(20, 3))
        path1 = Path(coordinates=coords)
        path2 = Path(coordinates=coords)
        path3 = path1 + path2
        assert path3 == path1 + path2
        assert np.allclose(
            path3.coordinates, np.concatenate((path1.coordinates, path2.coordinates))
        )


class TestRandomWalk(BaseTest):
    def test_extend_coordinates(self):
        path = Path()
        num_sites = NumSites(80)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path=path,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=14,
            chunk_size=50,
        )
        assert len(path.coordinates) == 80

    def test_default_args(self):
        path = hard_sphere_random_walk(termination=5)
        assert len(path.coordinates) == 5
        diffs = path.coordinates[0:-2] - path.coordinates[1:-1]
        assert np.allclose(0.15, np.linalg.norm(diffs, axis=1), atol=1e-4)

    def test_random_walk(self):
        path = Path()
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path=path,
            termination=[num_sites, max_attempts],
            bond_length=0.25,
            radius=0.22,
            seed=14,
        )
        assert len(path.coordinates) == 20
        diffs = path.coordinates[0:-2] - path.coordinates[1:-1]
        assert np.allclose(0.25, np.linalg.norm(diffs, axis=1), atol=1e-4)
        comp = path.to_compound()
        assert comp.n_particles == 20
        assert comp.n_bonds == 19
        # Test bounds of random initial point
        assert np.all(np.abs(path.coordinates[0])) < 20 * 0.22

    def test_dense_random_walk(self):
        rw_system = Path()
        vol_constaint = CuboidConstraint(center=(0, 0, 0), Lx=8)
        radius = 0.25
        bond_L = radius + 0.001
        num_chains = 200
        chain_lengths = 40
        tolerance = 0.01

        for i in range(num_chains):
            chain_passed = False
            attempt = 0
            while not chain_passed:
                initial_point = vol_constaint.find_low_density_points(
                    n_candidates=200 + attempt,
                    points=rw_system.coordinates,
                    buffer=radius,
                )
                term = Termination(
                    [NumSites(chain_lengths), NumAttempts(chain_lengths)]
                )
                try:
                    hard_sphere_random_walk(
                        path=rw_system,
                        bead_name="A",
                        radius=radius,
                        bond_length=bond_L,
                        volume_constraint=vol_constaint,
                        termination=term,
                        seed=i,
                        initial_point=initial_point[attempt],
                        rw_angles=AnglesSampler("normal", dict(loc=2.4, scale=1)),
                        tolerance=tolerance,
                    )
                    if term.success:
                        attempt = 0
                        chain_passed = True
                    else:
                        attempt += 1
                except PathConvergenceError:
                    attempt += 1  # try next initial_point candidate
                except Exception:
                    break  # only break on unexpected errors
        comp = rw_system.to_compound()
        assert len(comp.check_for_overlap(minimum_distance=radius - tolerance)) == 0

    def test_include_compound(self):
        L = 3
        box = mb.fill_box(
            compound=mb.load("C", smiles=True), n_compounds=20, box=[L, L, L]
        )
        box.translate_to((0, 0, 0))
        vol_constraint = CuboidConstraint(Lx=L, center=(0, 0, 0))
        path = Path()
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        chain = hard_sphere_random_walk(
            path=path,
            termination=[num_sites, max_attempts],
            bond_length=0.25,
            volume_constraint=vol_constraint,
            include_compound=box,
            radius=0.22,
            seed=14,
        )
        box.add(chain.to_compound())
        assert (
            len(box.check_for_overlap(excluded_bond_depth=2, minimum_distance=0.22))
            == 0
        )

    def test_set_initial_point(self):
        path = Path()
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path=path,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            initial_point=(1, 2, 3),
            seed=14,
        )
        assert np.array_equal(path.coordinates[0], np.array([1, 2, 3]))

    @pytest.mark.parametrize(
        "index", [2, np.int64(2), np.array(2), np.array([2])], ids=type
    )
    def test_initial_point_index_types(self, index):
        """Every spelling of a site index starts the same walk from that site."""

        def walk_from_index(initial_point):
            path = straight_line(spacing=0.25, N=5)
            hard_sphere_random_walk(
                path=path,
                termination=Termination([NumSites(10), NumAttempts(1e4)]),
                bond_length=0.25,
                radius=0.22,
                initial_point=initial_point,
                connectivity="link-linear",
                seed=14,
            )
            return path

        path = walk_from_index(index)
        # The first site of the walk is bonded to, and one bond length from, site 2.
        assert path.bond_graph.has_edge(2, 5)
        assert np.isclose(
            np.linalg.norm(path.coordinates[5] - path.coordinates[2]), 0.25
        )
        # Every index type gives an identical path.
        expected = walk_from_index(2)
        assert np.allclose(path.coordinates, expected.coordinates)
        assert set(path.bond_graph.edges) == set(expected.bond_graph.edges)

    def test_initial_point_bad_type(self):
        for bad in [2.5, np.array([1, 2]), "2"]:
            with pytest.raises(ValueError):
                hard_sphere_random_walk(
                    path=straight_line(spacing=0.25, N=5),
                    termination=Termination([NumSites(10), NumAttempts(1e4)]),
                    bond_length=0.25,
                    radius=0.22,
                    initial_point=bad,
                )

    def test_seeds(self):
        path1 = Path()
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path=path1,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=14,
        )
        path2 = Path()
        hard_sphere_random_walk(
            path=path2,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=14,
        )
        assert np.allclose(path1.coordinates, path2.coordinates, atol=1e-7)

    def test_seeds_user_angles_sampler(self):
        # A reused, stateful user sampler is still driven by the walk's seed,
        # so same-seed walks reproduce.
        sampler = AnglesSampler("normal", dict(loc=2.4, scale=0.3))
        kwargs = dict(bond_length=0.25, radius=0.22, rw_angles=sampler, seed=14)
        path1 = Path()
        hard_sphere_random_walk(
            path=path1,
            termination=Termination([NumSites(20), NumAttempts(1e4)]),
            **kwargs,
        )
        path2 = Path()
        hard_sphere_random_walk(
            path=path2,
            termination=Termination([NumSites(20), NumAttempts(1e4)]),
            **kwargs,
        )
        assert np.allclose(path1.coordinates, path2.coordinates, atol=1e-7)

    def test_seeds_random_namer(self):
        # An unseeded stochastic namer is driven by the walk's seed, so
        # same-seed walks reproduce both names and coordinates.
        kwargs = dict(radius=0.1, bond_length=0.15, termination=20, seed=42)
        path1 = hard_sphere_random_walk(bead_name=RandomNamer(["_A", "_B"]), **kwargs)
        path2 = hard_sphere_random_walk(bead_name=RandomNamer(["_A", "_B"]), **kwargs)
        assert list(path1.beads) == list(path2.beads)
        assert np.allclose(path1.coordinates, path2.coordinates, atol=1e-7)
        # Naming draws from a separate substream, so the namer choice does not
        # perturb geometry: same seed gives the same coordinates as a constant namer.
        const = hard_sphere_random_walk(bead_name="_A", **kwargs)
        assert np.allclose(path1.coordinates, const.coordinates, atol=1e-7)

    def test_from_path(self):
        path1 = Path()
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path=path1,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=24,
        )
        path2 = Path()
        hard_sphere_random_walk(
            path=path2,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=24,
        )
        assert len(path1.coordinates) == 20
        assert len(path2.coordinates) == 20
        assert np.allclose(path1.coordinates, path2.coordinates, atol=1e-6)

    @pytest.mark.parametrize(
        "constraint",
        [
            CuboidConstraint(Lx=6, Ly=6, Lz=6),
            SphereConstraint(center=(0, 0, 0), radius=3),
            CylinderConstraint(radius=3, height=6),
        ],
    )
    def test_seeds_with_volume_constraint(self, constraint):
        termination = Termination([NumSites(15), NumAttempts(1e4)])
        path1 = Path()
        hard_sphere_random_walk(
            path=path1,
            termination=termination,
            bond_length=0.25,
            radius=0.22,
            volume_constraint=constraint,
            seed=14,
        )
        path2 = Path()
        hard_sphere_random_walk(
            path=path2,
            termination=termination,
            bond_length=0.25,
            radius=0.22,
            volume_constraint=constraint,
            seed=14,
        )
        assert len(path1.coordinates) == len(path2.coordinates)
        assert np.allclose(path1.coordinates, path2.coordinates, atol=1e-6)

    def test_walk_inside_cube(self):
        path = Path()
        cube = CuboidConstraint(Lx=5, Ly=5, Lz=5)
        for i in range(100):
            hard_sphere_random_walk(
                path=path,
                termination=Termination([NumSites(5), NumAttempts(100)]),
                bond_length=0.25,
                radius=0.22,
                volume_constraint=cube,
                seed=14,
                tolerance=0.2,
            )
        bounds = bounding_box(path.coordinates)
        assert np.all(bounds < np.array([5 - 0.4, 5 - 0.4, 5 - 0.4]))

    def test_walk_inside_cube_with_pbc(self):
        # First make sure this seed gives a path outside these bounds without PBC
        path1 = Path()
        hard_sphere_random_walk(
            path=path1,
            termination=Termination([NumSites(500), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            initial_point=(0, 0, 0),
            volume_constraint=None,
            seed=14,
        )
        comp = path1.to_compound()
        assert np.all(comp.get_boundingbox().lengths > np.array([5, 5, 5]))

        path2 = Path()
        cube = CuboidConstraint(Lx=5, Ly=5, Lz=5, pbc=(True, True, True))
        hard_sphere_random_walk(
            path=path2,
            termination=Termination([NumSites(500), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            initial_point=(0, 0, 0),
            volume_constraint=cube,
            seed=14,
        )
        comp = path2.to_compound()
        assert np.all(comp.get_boundingbox().lengths <= np.array([5, 5, 5]))

    def test_walk_inside_sphere(self):
        path = Path()
        sphere = SphereConstraint(radius=4, center=(2, 2, 2))
        hard_sphere_random_walk(
            path=path,
            termination=Termination([NumSites(200), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            volume_constraint=sphere,
            initial_point=(0, 0, 0),
            seed=90,
        )
        bounds = bounding_box(path.coordinates)
        assert np.all(bounds < np.array([(2 * 4) - 0.22]))

    def test_walk_inside_cylinder(self):
        path = Path()
        cylinder = CylinderConstraint(radius=3, height=6, center=(0, 0, 0))
        hard_sphere_random_walk(
            path=path,
            termination=Termination([NumSites(200), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            volume_constraint=cylinder,
            seed=14,
        )
        bounds = bounding_box(path.coordinates)
        extents = bounds[1] - bounds[0]  # max - min for each dimension
        assert extents[0] < 6 - 0.22 * 2  # x extent
        assert extents[1] < 6 - 0.22 * 2  # y extent
        assert extents[2] < 6 - 0.22 * 2  # z extent

    @pytest.mark.parametrize("a_len, b_len", [(2, 3), (3, 5), (4, 4)])
    def test_multiple_paths(self, a_len, b_len):
        # Initialize two random walks
        chainDict = {
            "_A": {"bond_length": 0.15, "n_mers": a_len},
            "_B": {"bond_length": 0.15, "n_mers": b_len},
        }
        aPath = Path()

        # Randomize selection
        chainsList = ["_A", "_A", "_B", "_B", "_B"]

        # Fold together
        for chain in chainsList:
            # Create finish criteria
            num_sites = NumSites(chainDict[chain]["n_mers"])

            # Run a random walk
            hard_sphere_random_walk(
                path=aPath,
                radius=chainDict[chain]["bond_length"] * 0.95,
                bond_length=chainDict[chain]["bond_length"],
                termination=num_sites,
                bead_name=chain,
                seed=14,
            )

        data = aPath.bond_graph.nodes(data=True)
        site_idx = 0
        for chain in chainsList:
            assert data[site_idx]["name"] == chain
            site_idx += chainDict[chain]["n_mers"]

    def test_random_walk_bonds(self):
        path = hard_sphere_random_walk(
            radius=0.2, bond_length=0.25, termination=5, connectivity="linear"
        )
        assert len(path.bond_graph.edges) == 4
        hard_sphere_random_walk(
            path,
            radius=0.2,
            bond_length=0.25,
            termination=5,
            connectivity="link-linear",
        )
        assert len(path.bond_graph.edges) == 9
        # The new segment links to the last site of the previous chain.
        assert path.bond_graph.has_edge(4, 5)

        path = hard_sphere_random_walk(
            radius=0.2, bond_length=0.25, termination=5, connectivity="disconnected"
        )
        assert len(path.bond_graph.edges) == 0

        hard_sphere_random_walk(
            path, radius=0.2, bond_length=0.25, termination=5, connectivity="cycle"
        )
        assert len(path.bond_graph.edges) == 5

    def test_build_random_walk(self):
        path = Path()
        assert len(path.coordinates) == 0
        num_sites = NumSites(10)
        max_attempts = NumAttempts(1e4)

        hard_sphere_random_walk(
            path=path,
            radius=1.0,
            bond_length=1.1,
            termination=Termination([num_sites, max_attempts]),
        )
        assert len(path.coordinates) == 10
        # Check that bond graph has nodes
        assert path.bond_graph.number_of_nodes() == 10

    def test_init_positions(self):
        # Random point
        path = hard_sphere_random_walk(
            radius=1.0,
            bond_length=1.5,
            termination=1,
        )
        assert all([np.abs(coord) <= 1 for coord in path.coordinates[0]])

        # specific starting point
        path = hard_sphere_random_walk(
            radius=1.0,
            bond_length=1.5,
            initial_point=(0, 0, 0),
            termination=3,
        )
        assert np.allclose(path.coordinates[0], (0, 0, 0))
        assert len(path.coordinates) == 3

        # build from a previous path
        bond_length = 0.15
        hard_sphere_random_walk(
            path,
            radius=0.1,
            bond_length=bond_length,
            initial_point=0,
            termination=1,
        )
        assert np.allclose(np.linalg.norm(path.coordinates[3]), bond_length)
        assert len(path.coordinates) == 4

        # generate within a constraint
        constraint = CuboidConstraint.from_array([1, 1, 1])
        path = hard_sphere_random_walk(
            radius=0.1,
            bond_length=bond_length,
            volume_constraint=constraint,
            termination=1,
        )
        assert all(constraint.is_inside(points=path.coordinates, buffer=0.1))
        assert all([np.abs(coord) < 1 / 2 for coord in path.coordinates[0]])

        constraint = CuboidConstraint.from_array([1, 1, 1], center=(-0.5, -0.5, -0.5))
        path = hard_sphere_random_walk(
            radius=0.1,  # need a smaller buffer
            bond_length=0.2,
            volume_constraint=constraint,
            initial_point=(-0.25, -0.25, -0.25),
            termination=3,
            seed=100,
            tolerance=0.1,
        )
        assert np.allclose(path.coordinates[0], np.array([-0.25, -0.25, -0.25]))
        assert all(constraint.is_inside(points=path.coordinates[1:], buffer=0.1))
        for coord in path.coordinates[1:]:
            assert all([x < 0 for x in coord])

    def test_rw_normal_angles(self):
        from scipy.stats import normaltest

        num_sites = NumSites(1000)
        path = hard_sphere_random_walk(  # TODO: Map Angles into 0 to np.pi domain
            termination=num_sites,
            radius=0.0001,
            bond_length=1,
            rw_angles={
                "loc": np.pi / 2,
                "scale": 0.001,
            },  # larger scale doesn't center at mean
        )
        angles = []
        for i, j, k in zip(
            path.coordinates, path.coordinates[1:], path.coordinates[2:]
        ):
            BA = i - j
            BC = k - j
            norm_BA = np.linalg.norm(BA)
            norm_BC = np.linalg.norm(BC)
            angles.append(np.arccos(np.dot(BA, BC) / (norm_BA * norm_BC)))
        _, p_value = normaltest(angles)
        assert np.isclose(np.mean(angles), np.pi / 2, atol=1e-1)
        assert p_value > 0.05

    def test_rw_normal_angles_large_std(self):
        from scipy.stats import normaltest

        num_sites = NumSites(100)
        path = hard_sphere_random_walk(
            termination=num_sites,
            radius=0.001,  # point particle so radius doesn't influence selection
            bond_length=1,
            rw_angles={"loc": np.pi / 2, "scale": 0.5},
            trial_batch_size=8,
        )
        angles = []
        for i, j, k in zip(
            path.coordinates, path.coordinates[1:], path.coordinates[2:]
        ):
            BA = i - j
            BC = k - j
            norm_BA = np.linalg.norm(BA)
            norm_BC = np.linalg.norm(BC)
            angles.append(np.arccos(np.dot(BA, BC) / (norm_BA * norm_BC)))
        _, p_value = normaltest(angles)
        assert np.isclose(np.mean(angles), np.pi / 2, atol=1e-1)
        assert p_value > 0.05

    def test_rw_uniform_angles(self):
        import scipy.stats

        num_sites = NumSites(1000)
        min_max_angles = (np.pi / 3, np.pi / 2)
        path = hard_sphere_random_walk(
            termination=num_sites,
            radius=0.001,  # choose a point particle
            bond_length=1,
            rw_angles=min_max_angles,
            trial_batch_size=1,
        )
        angles = []
        for i, j, k in zip(
            path.coordinates, path.coordinates[1:], path.coordinates[2:]
        ):
            BA = i - j
            BC = k - j
            norm_BA = np.linalg.norm(BA)
            norm_BC = np.linalg.norm(BC)
            angles.append(np.arccos(np.dot(BA, BC) / (norm_BA * norm_BC)))
        uniform_loc_scale = (min_max_angles[0], min_max_angles[1] - min_max_angles[0])
        _, p_val = scipy.stats.kstest(angles, "uniform", args=uniform_loc_scale)
        assert p_val > 0.05
        assert np.isclose(np.mean(angles), np.pi * 5 / 12, atol=1e-1)

    def test_rw_cyclic_namer_sequence(self):
        path = hard_sphere_random_walk(
            radius=0.1,
            bond_length=0.15,
            termination=6,
            bead_name=CyclicNamer(["_A", "_B"]),
            seed=42,
        )
        assert list(path.beads) == ["_A", "_B", "_A", "_B", "_A", "_B"]

    def test_rw_cyclic_namer_blocks(self):
        path = hard_sphere_random_walk(
            radius=0.1,
            bond_length=0.15,
            termination=6,
            bead_name=CyclicNamer([("_A", 3), ("_B", 3)]),
            seed=42,
        )
        assert list(path.beads) == ["_A", "_A", "_A", "_B", "_B", "_B"]

    def test_rw_random_namer_draws_from_pool(self):
        path = hard_sphere_random_walk(
            radius=0.1,
            bond_length=0.15,
            termination=20,
            bead_name=RandomNamer(["_A", "_B"], seed=0),
            seed=42,
        )
        assert set(path.beads) <= {"_A", "_B"}
        assert len(path.beads) == 20

    def test_rw_namer_beads_in_bond_graph(self):
        path = hard_sphere_random_walk(
            radius=0.1,
            bond_length=0.15,
            termination=4,
            bead_name=CyclicNamer(["_A", "_B"]),
            seed=42,
        )
        node_names = [d["name"] for _, d in path.bond_graph.nodes(data=True)]
        assert node_names == ["_A", "_B", "_A", "_B"]

    def test_bond_length_less_than_radius(self):
        path = hard_sphere_random_walk(
            bead_name="A", bond_length=0.285, radius=0.392, termination=200, seed=7
        )
        coordinates = path.coordinates
        assert len(coordinates) == 200
        # Bonded neighbors sit at the bond length.
        bond_lengths = np.linalg.norm(np.diff(coordinates, axis=0), axis=1)
        assert np.allclose(bond_lengths, 0.285, atol=1e-5)
        # Every pair more than one bond apart respects the radius.
        distances = np.linalg.norm(
            coordinates[:, None, :] - coordinates[None, :, :], axis=-1
        )
        indices = np.arange(len(coordinates))
        non_bonded = np.abs(indices[:, None] - indices[None, :]) > 1
        assert distances[non_bonded].min() >= 0.392 - 1e-5

    @pytest.mark.parametrize(
        "bond_length, radius, rw_angles, expected",
        [
            (0.1, 0.25, None, "raise"),
            (0.15, 0.15, (0.3, 0.5), "raise"),
            (0.285, 0.392, (np.pi / 3, np.pi), "warn"),
        ],
    )
    def test_angle_range_against_radius(
        self, bond_length, radius, rw_angles, expected, caplog
    ):
        kwargs = dict(
            bond_length=bond_length,
            radius=radius,
            rw_angles=rw_angles,
            termination=30,
            seed=1,
        )
        if expected == "raise":
            with pytest.raises(ValueError):
                hard_sphere_random_walk(**kwargs)
        else:
            hard_sphere_random_walk(**kwargs)
            assert "overlap the site two bonds back" in caplog.text

    def test_unbonded_branch_start_uses_radius(self):
        # Under linear connectivity the first site of the walk is not bonded to
        # the site it starts from, so it is placed no closer than the radius.
        path = straight_line(spacing=0.5, N=6)
        hard_sphere_random_walk(
            path=path,
            bond_length=0.285,
            radius=0.392,
            termination=16,
            initial_point=2,
            connectivity="linear",
            seed=3,
        )
        assert not path.bond_graph.has_edge(2, 6)
        separation = np.linalg.norm(path.coordinates[6] - path.coordinates[2])
        assert separation >= 0.392 - 1e-5


class TestPathUtils(BaseTest):
    def test_target_sq_distances_no_pbc(self):
        target = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        points = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 3.0],
            ],
            dtype=np.float32,
        )
        d2 = target_sq_distances(target, points)
        expected = np.array([1.0, 4.0, 9.0], dtype=np.float32)
        assert np.allclose(d2, expected)

    def test_target_sq_distances_with_pbc(self):
        target = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        # These should wrap up -1
        points = np.array(
            [
                [9.0, 0.0, 0.0],
                [0.0, 9.0, 0.0],
                [0.0, 0.0, 9.0],
            ],
            dtype=np.float32,
        )
        pbc = np.array([True, True, True])
        box = np.array([10.0, 10.0, 10.0], dtype=np.float32)
        d2 = target_sq_distances(
            target,
            points,
            pbc=pbc,
            box_lengths=box,
        )
        expected = np.array([1.0, 1.0, 1.0], dtype=np.float32)
        assert np.allclose(d2, expected)

    def test_local_density_basic(self):
        candidate = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        targets = np.array(
            [
                [0.5, 0.0, 0.0],  # inside
                [1.5, 0.0, 0.0],  # outside
                [0.0, 0.5, 0.0],  # inside
            ],
            dtype=np.float32,
        )
        r_cut = 1.0
        density = local_density(candidate, targets, r_cut)
        assert density == 2

    def test_local_density_excludes_cutoff_boundary(self):
        candidate = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        targets = np.array(
            [
                [1.0, 0.0, 0.0],  # exactly r_cut
            ],
            dtype=np.float32,
        )
        r_cut = 1.0
        density = local_density(candidate, targets, r_cut)
        assert density == 0

    def test_local_density_no_targets(self):
        candidate = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        targets = np.empty((0, 3), dtype=np.float32)
        r_cut = 1.0
        density = local_density(candidate, targets, r_cut)
        assert density == 0

    def test_target_density_batch(self):
        candidates = np.array(
            [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ],
            dtype=np.float32,
        )
        targets = np.array(
            [
                [0.5, 0.0, 0.0],
                [1.5, 0.0, 0.0],
            ],
            dtype=np.float32,
        )
        r_cut = 1.0
        densities = target_density(candidates, targets, r_cut)
        # First candidate sees one neighbor
        # Second candidate sees one neighbor
        assert np.allclose(densities, np.array([1.0, 1.0], dtype=np.float32))

    def test_target_density_shape_and_dtype(self):
        candidates = np.random.rand(4, 3).astype(np.float32)
        targets = np.random.rand(10, 3).astype(np.float32)
        r_cut = 0.5
        densities = target_density(candidates, targets, r_cut)
        assert densities.shape == (4,)
        assert densities.dtype == np.float32

    @pytest.mark.parametrize(
        "distribution, kwargs, reference",
        [
            ("uniform", {"low": 0, "high": 1}, ("kstest", {"args": (0, 1)})),
            ("normal", {"loc": 0, "scale": 1}, ("normaltest", {})),
        ],
    )
    def test_angles_sampler(self, distribution, kwargs, reference):
        import scipy.stats

        from mbuild.path.points import AnglesSampler

        sampler = AnglesSampler(distribution, kwargs, rng=np.random.default_rng(0))
        points = sampler.sample(1000)
        if reference[0] == "kstest":
            _, p_value = scipy.stats.kstest(points, "uniform", **reference[1])
        elif reference[0] == "normaltest":
            _, p_value = scipy.stats.normaltest(points)

        assert p_value > 0.05

    @pytest.mark.parametrize("axis", [0, 1, 2])
    def test_check_path_pbc_overlap_across_face(self, axis):
        # Two beads sit just inside opposite faces of a 10 x 10 x 10 box.
        # minimum-image separation is 0.2
        box_lengths = np.full(3, 10.0, dtype=np.float32)
        radius = 1.0
        tolerance = 1e-5
        existing = np.full((1, 3), 5.0, dtype=np.float32)
        existing[0, axis] = 0.1
        candidate = np.full(3, 5.0, dtype=np.float32)
        candidate[axis] = 9.9
        # No PBCs should be accepted.
        assert check_path(existing, candidate, radius, tolerance)

        # Periodic on each axis, should be rejected.
        pbc = np.array([False, False, False])
        pbc[axis] = True
        assert not check_path(
            existing,
            candidate,
            radius,
            tolerance,
            pbc=pbc,
            box_lengths=box_lengths,
        )

    def test_check_path_excluded_indices(self):
        existing = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float32
        )
        candidate = np.array([1.1, 0.0, 0.0], dtype=np.float32)
        # Overlaps index 1 when every point is checked.
        assert not check_path(existing, candidate, 0.5, 1e-5)
        # Accepted once index 1 is excluded.
        assert check_path(
            existing,
            candidate,
            0.5,
            1e-5,
            excluded_indices=np.array([1], dtype=np.int64),
        )
        # Excluding a point that does not overlap changes nothing.
        assert not check_path(
            existing,
            candidate,
            0.5,
            1e-5,
            excluded_indices=np.array([0], dtype=np.int64),
        )
