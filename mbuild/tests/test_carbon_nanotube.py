import numpy as np
import pytest

from mbuild.lib.recipes import CarbonNanotube
from mbuild.tests.base_test import BaseTest


class TestCarbonNanotube(BaseTest):
    @pytest.mark.parametrize("n,m", [(5, 5), (10, 0), (6, 4)])
    def test_geometry(self, n, m):
        tube = CarbonNanotube(n=n, m=m, length=2.0)
        radii = np.linalg.norm(tube.xyz[:, :2], axis=1)
        assert np.allclose(radii, tube.radius)
        assert tube.xyz[:, 2].min() >= 0
        assert np.all(
            np.array(
                [len(list(tube.bond_graph.neighbors(p))) for p in tube.particles()]
            )
            <= 3
        )

    def test_periodic(self):
        tube = CarbonNanotube(n=5, m=5, length=3.0, periodic=True)
        assert tube.periodicity == (False, False, True)
        assert tube.n_bonds == 3 * tube.n_particles / 2

    def test_from_radius(self):
        tube = CarbonNanotube(radius=0.5, chirality="zigzag", length=1.0)
        assert tube.m == 0
        assert abs(tube.radius - 0.5) < 0.05

    def test_invalid_indices(self):
        with pytest.raises(ValueError):
            CarbonNanotube(n=3, m=5)
