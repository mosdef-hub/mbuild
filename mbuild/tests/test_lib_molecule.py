import numpy as np
import pytest

from mbuild.lib.molecules import (
    WaterSPC,
    WaterTIP3P,
    WaterTIP4P,
    WaterTIP4P2005,
    WaterTIP4PIce,
)
from mbuild.tests.base_test import BaseTest


class TestWater(BaseTest):
    @pytest.mark.parametrize(
        "model, bond_length, angle",
        [(WaterTIP3P, 0.09572, 104.52), (WaterSPC, 0.1, 109.47), (WaterOPC3, 0.09789, 109.47)],
    )
    def test_water_3site(self, model, bond_length, angle):
        water = model()
        o1 = [p for p in water.particles_by_name("OW")]
        assert len(o1) == 1
        o1 = o1[0]
        h1 = [p for p in water.particles_by_name("HW1")]
        assert len(h1) == 1
        h1 = h1[0]
        h2 = [p for p in water.particles_by_name("HW2")]
        assert len(h2) == 1
        h2 = h2[0]
        v1 = h1.xyz[0] - o1.xyz[0]
        v2 = h2.xyz[0] - o1.xyz[0]
        assert np.allclose(np.linalg.norm(v1), bond_length)
        assert np.allclose(np.linalg.norm(v2), bond_length)
        assert np.allclose(
            np.degrees(
                np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            ),
            angle,
        )

    @pytest.mark.parametrize(
        "model, bond_length, vsite_length, angle",
        [
            (WaterTIP4P, 0.09572, 0.015, 104.52),
            (WaterTIP4PIce, 0.09572, 0.01577, 104.52),
            (WaterTIP4P2005, 0.09572, 0.01546, 104.52),
            (WaterOPC, 0.08724, 0.01594, 103.6),
        ],
    )
    def test_water_4site(self, model, bond_length, vsite_length, angle):
        water = model()
        o1 = [p for p in water.particles_by_name("OW")]
        assert len(o1) == 1
        o1 = o1[0]
        h1 = [p for p in water.particles_by_name("HW1")]
        assert len(h1) == 1
        h1 = h1[0]
        h2 = [p for p in water.particles_by_name("HW2")]
        assert len(h2) == 1
        h2 = h2[0]
        m1 = [p for p in water.particles_by_name("MW")]
        assert len(m1) == 1
        m1 = m1[0]
        v1 = h1.xyz[0] - o1.xyz[0]
        v2 = h2.xyz[0] - o1.xyz[0]
        v3 = m1.xyz[0] - o1.xyz[0]
        assert np.allclose(np.linalg.norm(v1), bond_length)
        assert np.allclose(np.linalg.norm(v2), bond_length)
        assert np.allclose(np.linalg.norm(v3), vsite_length)
        assert np.allclose(
            np.degrees(
                np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            ),
            angle,
        )
