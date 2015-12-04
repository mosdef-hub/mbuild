import mbuild as mb
from mbuild.tests.base_test import BaseTest

from mbuild.lib.surfaces import Betacristobalite
from mbuild.lib.atoms import H


class TestMonolayer(BaseTest):

    def test_monolayer(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)

        chain = mb.Polymer(ch2, n=10)
        monolayer = mb.Monolayer(surface=Betacristobalite(), chain=chain,
                                 backfill=H(), pattern=pattern)

        assert monolayer.n_particles == 1800 + n * m * (10*3) + (100 - n*m)
        assert monolayer.n_bonds == 2300 + n * m * (10 * 2 + 9 + 1) + (100 - n * m)

    def test_pattern_kwargs(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)

        chain = mb.Polymer(ch2, n=10)
        monolayer = mb.Monolayer(surface=Betacristobalite(), chain=H(),
                                 guest_port_name='up', backfill=chain,
                                 backfill_port_name='down', pattern=pattern)

        chains = 100 - (n*m)

        assert monolayer.n_particles == 1800 + chains * (10*3) + (100 - chains)
        assert monolayer.n_bonds == 2300 + chains * (10 * 2 + 9 + 1) + (100 - chains)
