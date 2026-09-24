import numpy as np
import pytest

import mbuild as mb
from mbuild.lib.recipes import CarbonNanotube
from mbuild.tests.base_test import BaseTest


class TestCarbonNanotube(BaseTest):
    @pytest.mark.parametrize("n,m", [(5, 5), (10, 0), (6, 4)])
    def test_geometry(self, n, m):
        tube = CarbonNanotube(n=n, m=m, length=2.0, cap=None)
        radii = np.linalg.norm(tube.xyz[:, :2], axis=1)
        assert np.allclose(radii, tube.radius)
        assert tube.xyz[:, 2].min() >= 0
        assert np.all(
            np.array(
                [len(list(tube.bond_graph.neighbors(p))) for p in tube.particles()]
            )
            <= 3
        )

    @pytest.mark.parametrize("n,m", [(5, 5), (10, 0), (6, 4)])
    def test_capped_ends(self, n, m):
        tube = CarbonNanotube(n=n, m=m, length=2.0, cap="F")
        assert not tube.cap_ports
        assert not tube.all_ports()
        for particle in tube.particles():
            expected = 3 if particle.name == "C" else 1
            assert len(list(tube.bond_graph.neighbors(particle))) == expected
        fluorines = list(tube.particles_by_name("F"))
        assert fluorines
        for f in fluorines:
            carbon = next(iter(tube.bond_graph.neighbors(f)))
            assert np.isclose(np.linalg.norm(f.pos - carbon.pos), 0.135)

    def test_uncapped_ports(self):
        tube = CarbonNanotube(n=6, m=6, length=1.0, cap=None)
        edge_carbons = [
            p for p in tube.particles() if len(list(tube.bond_graph.neighbors(p))) == 2
        ]
        assert len(tube.cap_ports) == len(edge_carbons)
        assert len(tube.all_ports()) == len(edge_carbons)

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

    def test_invalid_cap(self):
        with pytest.raises(ValueError):
            CarbonNanotube(n=5, m=5, cap="Xe")

    def test_compound_cap(self):
        hydroxyl = mb.load("O", smiles=True)
        hydroxyl.remove(hydroxyl.children[-1])
        hydroxyl.add(hydroxyl.all_ports()[0], label="up", containment=False)

        uncapped = CarbonNanotube(n=6, m=0, length=3.0, periodic=False, cap=None)
        n_ports = len(uncapped.cap_ports)

        tube = CarbonNanotube(n=6, m=0, length=3.0, periodic=False, cap=hydroxyl)
        assert not tube.cap_ports
        assert not tube.all_ports()
        oxygens = list(tube.particles_by_name("O"))
        assert len(oxygens) == n_ports
        assert len(list(tube.particles_by_name("H"))) == n_ports
        for particle in tube.particles():
            expected = {"C": 3, "O": 2, "H": 1}[particle.name]
            assert len(list(tube.bond_graph.neighbors(particle))) == expected
        for oxygen in oxygens:
            carbon = next(n for n in tube.bond_graph.neighbors(oxygen) if n.name == "C")
            # The C-O distance is set by the two port separations, not a lookup
            assert 0.10 < np.linalg.norm(oxygen.pos - carbon.pos) < 0.16
