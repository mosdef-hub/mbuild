import networkx as nx
import numpy as np
import pytest

import mbuild as mb
from mbuild.path import (
    Cyclic,
    HardSphereRandomWalk,
    Knot,
    Lamellar,
    Path,
    Spiral2D,
    StraightLine,
)
from mbuild.path.bias import TargetCoordinate
from mbuild.path.path_utils import (
    local_density,
    target_density,
    target_sq_distances,
)
from mbuild.path.termination import (
    EndToEndDistance,
    NumAttempts,
    NumSites,
    RadiusOfGyration,
    Termination,
    WithinCoordinate,
)
from mbuild.tests.base_test import BaseTest, radius_of_gyration
from mbuild.utils.geometry import bounding_box
from mbuild.utils.volumes import (
    CuboidConstraint,
    CylinderConstraint,
    SphereConstraint,
)


class TestPaths(BaseTest):
    def test_from_coordinates(self):
        coords = np.random.uniform(-5, 5, size=(20, 3))
        path = Path.from_coordinates(coordinates=coords, bond_graph=nx.Graph())
        assert np.array_equal(coords, path.coordinates)
        for idx, node in enumerate(path.bond_graph.nodes(data=True)):
            assert np.array_equal(node[1]["xyz"], coords[idx])

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
        path = StraightLine(spacing=0.20, N=5, direction=(1, 0, 0))
        assert len(path.coordinates) == 5
        assert path.bond_graph.number_of_edges() == 4
        # 5 sites = 4 bonds at 0.20 each
        assert np.linalg.norm(path.coordinates[-1] - path.coordinates[0]) == 0.80
        for edge in path.bond_graph.edges(data=True):
            assert np.array_equal(edge[2]["direction"], np.array([1, 0, 0]))

    def test_cyclic_parameters(self):
        path = Cyclic(spacing=1, N=20)
        # C = 2*pi*r
        radius = 10 / np.pi
        assert np.allclose(path.radius, radius, 1e-2)

        path = Cyclic(spacing=1, radius=10 / np.pi, N=None)
        assert path.N == 20

        path = Cyclic(N=20, radius=10 / np.pi)
        assert path.spacing == 1

    def test_cyclic_bonding(self):
        path = Cyclic(spacing=1, N=20)
        assert path.bond_graph.number_of_edges() == 20
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles

    def test_knot(self):
        path = Knot(spacing=0.25, N=50, m=3)
        assert path.bond_graph.number_of_edges() == 50
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles

    def test_knot_bad_arg(self):
        with pytest.raises(ValueError):
            Knot(spacing=0.25, N=50, m=2)

    def test_spiral(self):
        path = Spiral2D(spacing=0.25, N=50, a=0.5, b=2)
        assert path.bond_graph.number_of_edges() == 49
        comp = path.to_compound()
        assert comp.n_bonds == comp.n_particles - 1

    def test_lamellar(self):
        path = Lamellar(
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
        assert Lx == Lz  # stacking and layering directions
        assert Ly > Lx


class TestRandomWalk(BaseTest):
    def test_no_termination(self):
        with pytest.raises(RuntimeError):
            HardSphereRandomWalk(
                termination=None,
                bond_length=0.25,
                radius=0.22,
                min_angle=np.pi / 4,
                max_angle=np.pi,
                max_attempts=1e4,
                seed=14,
            )

    def test_rg_termination(self):
        num_sites = NumSites(20)
        rg = RadiusOfGyration(3)
        max_attempts = NumAttempts(1e4)
        termination = Termination([num_sites, rg, max_attempts])
        rw_path = HardSphereRandomWalk(
            termination=termination,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        assert np.allclose(radius_of_gyration(rw_path.coordinates), 3, atol=1e-1)

    def test_re_termination(self):
        num_sites = NumSites(20)
        re = EndToEndDistance(2)
        max_attempts = NumAttempts(1e4)
        termination = Termination([num_sites, re, max_attempts])
        rw_path = HardSphereRandomWalk(
            termination=termination,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        dist = np.linalg.norm(rw_path.coordinates[-1] - rw_path.coordinates[0])
        assert np.allclose(dist, 2, atol=1e-1)

    def test_within_coord_termination(self):
        bias = TargetCoordinate(target_coordinate=(2, 2, 2), weight=0.8)
        termination = Termination(
            [
                WithinCoordinate(target_coordinate=(2, 2, 2), distance=0.22),
                NumAttempts(1e3),
            ]
        )
        rw_path = HardSphereRandomWalk(
            termination=termination,
            bias=bias,
            initial_point=(0, 0, 0),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            trial_batch_size=100,
            seed=14,
        )
        dist = np.linalg.norm(rw_path.coordinates[-1] - np.array([2, 2, 2]))
        assert dist <= 0.22

        termination = Termination(
            [
                WithinCoordinate(
                    target_coordinate=(3, 3, 3), distance=0.0, tolerance=1e-1
                ),
                NumAttempts(100),
            ]
        )
        bias = TargetCoordinate(target_coordinate=(3, 3, 3), weight=1.0)
        rw_path = HardSphereRandomWalk(
            termination=termination,
            bias=bias,
            initial_point=(0, 0, 0),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=5e4,
            trial_batch_size=200,
            seed=14,
        )
        dist = np.linalg.norm(rw_path.coordinates[-1] - np.array([3, 3, 3]))
        assert np.allclose(dist, 0, atol=1e-1)

    def test_extend_coordinates(self):
        num_sites = NumSites(80)
        max_attempts = NumAttempts(1e4)
        rw_path = HardSphereRandomWalk(
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
            chunk_size=50,
        )
        assert len(rw_path.coordinates) == 80

    def test_random_walk(self):
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        rw_path = HardSphereRandomWalk(
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        assert len(rw_path.coordinates) == 20
        diffs = rw_path.coordinates[0:-2] - rw_path.coordinates[1:-1]
        assert np.allclose(0.25, np.linalg.norm(diffs, axis=1), atol=1e-4)
        comp = rw_path.to_compound()
        assert comp.n_particles == 20
        assert comp.n_bonds == 19
        # Test bounds of random initial point
        assert np.all(rw_path.coordinates[0]) < 20 * 0.22

    def test_set_initial_point(self):
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        rw_path = HardSphereRandomWalk(
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            initial_point=(1, 2, 3),
            seed=14,
        )
        assert np.array_equal(rw_path.coordinates[0], np.array([1, 2, 3]))

    def test_seeds(self):
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        rw_path_1 = HardSphereRandomWalk(
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        rw_path_2 = HardSphereRandomWalk(
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        assert np.allclose(rw_path_1.coordinates, rw_path_2.coordinates, atol=1e-7)

    def test_from_path(self):
        num_sites = NumSites(20)
        max_attempts = NumAttempts(1e4)
        rw_path = HardSphereRandomWalk(
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        rw_path2 = HardSphereRandomWalk(
            termination=Termination([num_sites, max_attempts]),
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=24,
            start_from_path=rw_path,
            start_from_path_index=-1,
        )
        assert len(rw_path.coordinates) == 20
        assert len(rw_path2.coordinates) == 40
        for coord1, coord2 in zip(rw_path.coordinates[:10], rw_path2.coordinates[:10]):
            assert np.allclose(coord1, coord2, atol=1e-6)

    def test_walk_inside_cube(self):
        cube = CuboidConstraint(Lx=5, Ly=5, Lz=5)
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(500), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            volume_constraint=cube,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        bounds = bounding_box(rw_path.coordinates)
        assert np.all(bounds < np.array([5 - 0.44, 5 - 0.44, 5 - 0.44]))

    def test_walk_inside_cube_with_pbc(self):
        # First make sure this seed gives a path outside these bounds without PBC
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(500), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            initial_point=(0, 0, 0),
            volume_constraint=None,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        comp = rw_path.to_compound()
        assert np.all(comp.get_boundingbox().lengths > np.array([5, 5, 5]))

        cube = CuboidConstraint(Lx=5, Ly=5, Lz=5, pbc=(True, True, True))
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(500), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            initial_point=(0, 0, 0),
            volume_constraint=cube,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        comp = rw_path.to_compound()
        assert np.all(comp.get_boundingbox().lengths <= np.array([5, 5, 5]))

    def test_walk_inside_sphere(self):
        sphere = SphereConstraint(radius=4, center=(2, 2, 2))
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(200), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            volume_constraint=sphere,
            initial_point=(0, 0, 0),
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=90,
        )
        bounds = bounding_box(rw_path.coordinates)
        assert np.all(bounds < np.array([(2 * 4) - 0.22]))

    def test_walk_inside_cylinder(self):
        cylinder = CylinderConstraint(radius=3, height=6, center=(0, 0, 0))
        rw_path = HardSphereRandomWalk(
            termination=Termination([NumSites(200), NumAttempts(1e4)]),
            bond_length=0.25,
            radius=0.22,
            volume_constraint=cylinder,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        bounds = bounding_box(rw_path.coordinates)
        assert bounds[0][0] < 6 - 0.22 * 2
        assert bounds[1][1] < 6 - 0.22 * 2
        assert bounds[2][2] < 6 - 0.22 * 2


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
