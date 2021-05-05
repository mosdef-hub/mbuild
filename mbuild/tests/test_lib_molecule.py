import numpy as np

from mbuild.tests.base_test import BaseTest
from mbuild.lib.molecules import WaterTIP3P
from mbuild.lib.molecules import WaterTIP4P
from mbuild.lib.molecules import WaterTIP4PIce
from mbuild.lib.molecules import WaterTIP4P2005
from mbuild.lib.molecules import WaterSPC


class TestWater(BaseTest):

    def test_tip3p(self):
        water = WaterTIP3P()
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
        assert np.allclose(np.linalg.norm(v1), 0.09572)
        assert np.allclose(np.linalg.norm(v2), 0.09572)
        assert np.allclose(
            np.degrees(
                np.arccos(
                    np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
                )
            ), 104.52
        )

    def test_spc(self):
        water = WaterSPC()
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
        assert np.allclose(np.linalg.norm(v1), 0.1)
        assert np.allclose(np.linalg.norm(v2), 0.1)
        assert np.allclose(
            np.degrees(
                np.arccos(
                    np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
                )
            ), 109.47
        )

    def test_tip4p(self):
        water = WaterTIP4P()
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
        assert np.allclose(np.linalg.norm(v1), 0.09572)
        assert np.allclose(np.linalg.norm(v2), 0.09572)
        assert np.allclose(np.linalg.norm(v3), 0.015)
        assert np.allclose(
            np.degrees(
                np.arccos(
                    np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
                )
            ), 104.52
        )

    def test_tip4pice(self):
        water = WaterTIP4PIce()
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
        assert np.allclose(np.linalg.norm(v1), 0.09572)
        assert np.allclose(np.linalg.norm(v2), 0.09572)
        assert np.allclose(np.linalg.norm(v3), 0.01577)
        assert np.allclose(
            np.degrees(
                np.arccos(
                    np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
                )
            ), 104.52
        )

    def test_tip4p2005(self):
        water = WaterTIP4P2005()
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
        assert np.allclose(np.linalg.norm(v1), 0.09572)
        assert np.allclose(np.linalg.norm(v2), 0.09572)
        assert np.allclose(np.linalg.norm(v3), 0.01546)
        assert np.allclose(
            np.degrees(
                np.arccos(
                    np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
                )
            ), 104.52
        )
