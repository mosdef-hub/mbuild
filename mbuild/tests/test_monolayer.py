import mbuild as mb
from mbuild.lib.atoms import H
from mbuild.lib.recipes import Monolayer, Polymer
from mbuild.lib.surfaces import Betacristobalite
from mbuild.tests.base_test import BaseTest


class TestMonolayer(BaseTest):
    def test_monolayer(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)

        chain = Polymer(ch2, n=10)
        bc = Betacristobalite()
        monolayer = Monolayer(
            surface=bc, chains=chain, backfill=H(), pattern=pattern
        )

        assert monolayer.n_particles == 2000 + n * m * 29
        assert monolayer.n_bonds == 2500 + n * m * 29

    def test_pattern_kwargs(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)

        chain = Polymer(ch2, n=10)
        bc = Betacristobalite()
        monolayer = Monolayer(
            surface=bc,
            chains=H(),
            guest_port_name="up",
            backfill=chain,
            backfill_port_name="down",
            pattern=pattern,
        )

        chains = 100 - (n * m)

        assert monolayer.n_particles == 2000 + chains * 29
        assert monolayer.n_bonds == 2500 + chains * 29

    def test_mixed_monolayer(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)
        fractions = [0.75, 0.25]

        chain_a = Polymer(ch2, n=5)
        chain_b = Polymer(ch2, n=15)
        bc = Betacristobalite()
        monolayer = Monolayer(
            surface=bc,
            chains=[chain_a, chain_b],
            fractions=fractions,
            backfill=H(),
            pattern=pattern,
        )

        n_a = round(n * m * 0.75)
        n_b = round(n * m * 0.25)
        assert monolayer.n_particles == 2000 + n_a * 14 + n_b * 44
        assert monolayer.n_bonds == 2500 + n_a * 14 + n_b * 44
