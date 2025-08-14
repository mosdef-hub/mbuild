import networkx as nx
import numpy as np

from mbuild.path import HardSphereRandomWalk
from mbuild.tests.base_test import BaseTest
from mbuild.utils.geometry import bounding_box
from mbuild.utils.volumes import (
    CuboidConstraint,
    CylinderConstraint,
    SphereConstraint,
)


class TestRandomWalk(BaseTest):
    def test_random_walk(self):
        rw_path = HardSphereRandomWalk(
            N=20,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
            bond_graph=nx.path_graph(20),
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
        rw_path = HardSphereRandomWalk(
            N=20,
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
        rw_path_1 = HardSphereRandomWalk(
            N=20,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        rw_path_2 = HardSphereRandomWalk(
            N=20,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        assert np.allclose(rw_path_1.coordinates, rw_path_2.coordinates, atol=1e-7)

    def test_from_path(self):
        rw_path = HardSphereRandomWalk(
            N=20,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        rw_path2 = HardSphereRandomWalk(
            N=20,
            bond_length=0.25,
            radius=0.22,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
            start_from_path=rw_path,
            start_from_path_index=-1,
        )
        assert len(rw_path2.coordinates) == 40
        assert np.array_equal(rw_path.coordinates, rw_path2.coordinates[:20])

    def test_walk_inside_cube(self):
        cube = CuboidConstraint(Lx=5, Ly=5, Lz=5)
        rw_path = HardSphereRandomWalk(
            N=500,
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

    def test_walk_inside_sphere(self):
        sphere = SphereConstraint(radius=4, center=(2, 2, 2))
        rw_path = HardSphereRandomWalk(
            N=200,
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
        cylinder = CylinderConstraint(radius=3, height=6, center=(1.5, 1.5, 3))
        rw_path = HardSphereRandomWalk(
            N=200,
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
