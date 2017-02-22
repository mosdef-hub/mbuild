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
        monolayer = mb.Monolayer(surface=Betacristobalite(), chains=chain,
                                 backfill=H(), pattern=pattern)

        assert monolayer.n_particles == 1900 + n * m * (10*3) + (100 - n*m)
        assert monolayer.n_bonds == 2400 + n * m * (10 * 2 + 9 + 1) + (100 - n * m)

    def test_pattern_kwargs(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)

        chain = mb.Polymer(ch2, n=10)
        monolayer = mb.Monolayer(surface=Betacristobalite(), chains=H(),
                                 guest_port_name='up', backfill=chain,
                                 backfill_port_name='down', pattern=pattern)

        chains = 100 - (n*m)

        assert monolayer.n_particles == 1900 + chains * (10*3) + (100 - chains)
        assert monolayer.n_bonds == 2400 + chains * (10 * 2 + 9 + 1) + (100 - chains)

    def test_mixed_monolayer(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)
        fractions = [0.75,0.25]

        chain_a = mb.Polymer(ch2, n=5)
        chain_b = mb.Polymer(ch2, n=15)
        monolayer = mb.Monolayer(surface=Betacristobalite(),
                                 chains=[chain_a, chain_b],
                                 fractions=fractions,
                                 backfill=H(),
                                 pattern=pattern)

        n_a = round(n * m * 0.75)
        n_b = round(n * m * 0.25)
        assert monolayer.n_particles == 1900 + n_a * 5 * 3 + n_b * 15 * 3 + (100 - (n_a + n_b))
        assert monolayer.n_bonds == 2400 + n_a * (5 * 2 + 4 + 1) + n_b * (15 * 2 + 14 + 1) + (100 - (n_a + n_b))
