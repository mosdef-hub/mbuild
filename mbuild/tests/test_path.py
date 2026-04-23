import numpy as np
import pytest

import mbuild as mb
from mbuild.path.build import (
    Path,
    crosslink,
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
from mbuild.path.path_utils import (
    local_density,
    target_density,
    target_sq_distances,
)
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
        straight_line(path, spacing=0.20, N=5, direction=(1, 0, 0))
        assert len(path.coordinates) == 5
        assert path.bond_graph.number_of_edges() == 4
        # 5 sites = 4 bonds at 0.20 each
        assert np.allclose(
            np.linalg.norm(path.coordinates[-1] - path.coordinates[0]), 0.80
        )
        for edge in path.bond_graph.edges(data=True):
            assert np.allclose(edge[2]["direction"], np.array([1, 0, 0]))

    def test_cyclic_parameters(self):
        path = Path()
        cyclic(path, spacing=1, N=20)
        # C = 2*pi*r
        radius = 20 * 1 / (2 * np.pi)
        # Check that computed radius matches expected
        actual_radius = np.linalg.norm(path.coordinates[0])
        assert np.allclose(actual_radius, radius, atol=1e-2)

        path2 = Path()
        cyclic(path2, spacing=1, radius=10 / np.pi, N=None)
        assert len(path2.coordinates) == 20

        path3 = Path()
        cyclic(path3, N=20, radius=10 / np.pi, spacing=None)
        # Check spacing
        dist = np.linalg.norm(path3.coordinates[1] - path3.coordinates[0])
        assert np.allclose(dist, 1.0, atol=1e-2)

    def test_cyclic_bonding(self):
        path = Path()
        cyclic(path, spacing=1, N=20)
        assert path.bond_graph.number_of_edges() == 20
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles

    def test_knot(self):
        path = Path()
        knot(path, spacing=0.25, N=50, m=3)
        assert path.bond_graph.number_of_edges() == 50
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles

    def test_knot_bad_arg(self):
        path = Path()
        with pytest.raises(ValueError):
            knot(path, spacing=0.25, N=50, m=2)

    def test_spiral(self):
        path = Path()
        spiral_2D(path, N=50, a=0.5, b=2, spacing=0.25)
        assert path.bond_graph.number_of_edges() == 49
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles - 1

    def test_lamellar(self):
        path = Path()
        lamellar(
            path,
            bond_length=0.25,
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
            path_left_to_right,
            bond_length=0.25,
            num_layers=3,
            layer_separation=1.0,
            layer_length=3.0,
            num_stacks=3,
            stack_separation=1.0,
            initial_point=(0, 0, 0),
        )

        path_right_to_left = Path()
        lamellar(
            path_right_to_left,
            bond_length=0.25,
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
            path,
            bond_length=0.25,
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
        helix(path, N=50, radius=2.0, rise=0.5, twist=30)
        assert len(path.coordinates) == 50
        assert path.bond_graph.number_of_edges() == 49
        # Check that all points are roughly at radius distance from z-axis
        radii = np.sqrt(path.coordinates[:, 0] ** 2 + path.coordinates[:, 1] ** 2)
        assert np.allclose(radii, 2.0, atol=1e-5)

    def test_zigzag(self):
        path = Path()
        zigzag(path, N=20, spacing=1.0, angle_deg=120.0, sites_per_segment=5)
        assert len(path.coordinates) == 20
        assert path.bond_graph.number_of_edges() == 19


class TestRandomWalk(BaseTest):
    def test_extend_coordinates(self):
        path = Path()
        num_sites = NumSites(80)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path,
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
            path,
            termination=[num_sites, max_attempts],
            bond_length=0.25,
            radius=0.22,
            seed=14,
            run_on_gpu=True,
        )
        assert len(path.coordinates) == 20
        diffs = path.coordinates[0:-2] - path.coordinates[1:-1]
        assert np.allclose(0.25, np.linalg.norm(diffs, axis=1), atol=1e-4)
        comp = path.to_compound()
        assert comp.n_particles == 20
        assert comp.n_bonds == 19
        # Test bounds of random initial point
        assert np.all(np.abs(path.coordinates[0])) < 20 * 0.22

    def test_set_initial_point(self):
        path = Path()
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            initial_point=(1, 2, 3),
            seed=14,
        )
        assert np.array_equal(path.coordinates[0], np.array([1, 2, 3]))

    def test_seeds(self):
        path1 = Path()
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path1,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=14,
        )
        path2 = Path()
        hard_sphere_random_walk(
            path2,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=14,
        )
        assert np.allclose(path1.coordinates, path2.coordinates, atol=1e-7)

    def test_from_path(self):
        path1 = Path()
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        hard_sphere_random_walk(
            path1,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=24,
        )
        path2 = Path()
        hard_sphere_random_walk(
            path2,
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            seed=24,
        )
        assert len(path1.coordinates) == 20
        assert len(path2.coordinates) == 20
        assert np.allclose(path1.coordinates, path2.coordinates, atol=1e-6)

    def test_walk_inside_cube(self):
        path = Path()
        cube = CuboidConstraint(Lx=5, Ly=5, Lz=5)
        hard_sphere_random_walk(
            path,
            termination=Termination([NumSites(500), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            volume_constraint=cube,
            seed=14,
        )
        bounds = bounding_box(path.coordinates)
        assert np.all(bounds < np.array([5 - 0.44, 5 - 0.44, 5 - 0.44]))

    def test_walk_inside_cube_with_pbc(self):
        # First make sure this seed gives a path outside these bounds without PBC
        path1 = Path()
        hard_sphere_random_walk(
            path1,
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
            path2,
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
            path,
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
            path,
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
        assert len(path.coordinates == 3)

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
        assert len(path.coordinates == 4)

        # generate within a constraint
        constraint = CuboidConstraint.from_array([1, 1, 1])
        path = hard_sphere_random_walk(
            radius=0.1,
            bond_length=bond_length,
            volume_constraint=constraint,
            termination=1,
        )
        assert constraint.is_inside(points=path.coordinates, buffer=0.1)
        assert all([np.abs(coord) < 1 / 2 for coord in path.coordinates[0]])

        constraint = CuboidConstraint.from_array([1, 1, 1], center=(-0.5, -0.5, -0.5))
        path = hard_sphere_random_walk(
            radius=0.1,  # need a smaller buffer
            bond_length=0.2,
            volume_constraint=constraint,
            initial_point=(0, 0, 0),
            termination=3,
            seed=100,
        )
        assert np.allclose(path.coordinates[0], np.array([0, 0, 0]))
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

        sampler = AnglesSampler(distribution, kwargs, seed=0)
        points = sampler.sample(1000)
        if reference[0] == "kstest":
            _, p_value = scipy.stats.kstest(points, "uniform", **reference[1])
        elif reference[0] == "normaltest":
            _, p_value = scipy.stats.normaltest(points)

        assert p_value > 0.05


class TestCrossLinks(BaseTest):
    def test_find_links_line(self):
        path = Path()
        pos1 = np.zeros((10, 3))
        pos1[:, 1] = np.arange(10)
        path.append_coordinates(pos1)
        path.form_linear_bond_graph()

        pos2 = np.zeros((10, 3))
        pos2[:, 0] += 1
        pos2[:, 1] = np.arange(10)
        path.append_coordinates(pos2)
        path.form_linear_bond_graph(indices=np.arange(10, 20))

        for i in range(10):
            crosslink(
                path, initial_point=i, radius=1.1, excluded_bond_depth=10, seed=42
            )

        clinks = sum([b == "_R" for b in path.beads])
        bbones = sum([b == "_A" for b in path.beads])
        assert clinks == 10
        assert bbones == 20
        for i in range(10):
            assert (i, i + 20) in path.bond_graph.edges
            assert (i + 10, i + 20) in path.bond_graph.edges

    def test_deterministic_rw(self):
        path1 = hard_sphere_random_walk(  # TODO: hsrw cuts off some chains early
            radius=1,
            bond_length=2,
            termination=20,
            rw_angles=(np.pi / 2, np.pi),
            seed=1,
        )
        path2 = hard_sphere_random_walk(  # TODO: hsrw cuts off some chains early
            radius=1,
            bond_length=2,
            termination=20,
            rw_angles=(np.pi / 2, np.pi),
            seed=1,
        )
        assert path1 == path2

    def test_deterministic_crosslink(self):
        """Test set coordinates, relax seed as well"""
        rng = np.random.default_rng(1)

        points = rng.random((100, 3))
        path1 = Path(coordinates=points)
        path2 = Path(coordinates=points)
        assert path1 == path2

        crosslink(path1, radius=1, excluded_bond_depth=2)
        crosslink(path2, radius=1, excluded_bond_depth=2)
        assert path1 == path2

        path1.relax(0.2, None, steps=10)
        path2.relax(0.2, None, steps=10)
        print(path1.coordinates - path2.coordinates)
        assert np.allclose(path1.coordinates, path2.coordinates, atol=1e-6)
