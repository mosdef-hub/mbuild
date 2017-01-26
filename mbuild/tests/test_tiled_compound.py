import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestTiledCompound(BaseTest):

    def test_2d_replication(self, betacristobalite):
        nx = 2
        ny = 3
        nz = 1
        tiled = mb.TiledCompound(betacristobalite, [nx, ny, nz])
        assert tiled.n_particles == 1900 * nx * ny
        assert tiled.n_bonds == 2400 * nx * ny
        for at in tiled.particles():
            if at.name.startswith('Si'):
                assert len(tiled.bond_graph.neighbors(at)) <= 4
            elif at.name.startswith('O'):
                assert len(tiled.bond_graph.neighbors(at)) <= 2

        for at in tiled.particles():
            if at.name.startswith('Si'):
                assert len(tiled.bond_graph.neighbors(at)) <= 4
            elif at.name.startswith('O'):
                assert len(tiled.bond_graph.neighbors(at)) <= 2

    def test_no_replication(self, betacristobalite):
        nx = 1
        ny = 1
        nz = 1
        tiled = mb.TiledCompound(betacristobalite, [nx, ny, nz])
        assert tiled.n_particles == 1900 * nx * ny
        assert tiled.n_bonds == 2400 * nx * ny

    def test_incorrect_periodicity(self, betacristobalite):
        nx = 2
        ny = 3
        nz = 2
        with pytest.raises(ValueError):
            mb.TiledCompound(betacristobalite, [nx, ny, nz])

    def test_negative_periodicity(self, betacristobalite):
        nx = -2
        ny = 3
        nz = 2
        with pytest.raises(ValueError):
            mb.TiledCompound(betacristobalite, [nx, ny, nz])


