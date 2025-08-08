import numpy as np

from mbuild.path import HardSphereRandomWalk
from mbuild.tests.base_test import BaseTest
from mbuild.utils.geometry import bounding_box
from mbuild.utils.volumes import CuboidConstraint


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
        )
        assert len(rw_path.coordinates) == 20

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
        assert np.allclose(rw_path_1.coordinates, rw_path_2.coordinates, atol=1e-4)

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

    def test_walk_inside_cube(self):
        cube = CuboidConstraint(Lx=5, Ly=5, Lz=5)
        rw_path = HardSphereRandomWalk(
            N=50,
            bond_length=0.25,
            radius=0.22,
            volume_constraint=cube,
            min_angle=np.pi / 4,
            max_angle=np.pi,
            max_attempts=1e4,
            seed=14,
        )
        bounds = bounding_box(rw_path.coordinates)
        assert np.all(bounds < np.array([5, 5, 5]))
