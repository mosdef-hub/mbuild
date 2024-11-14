import pytest

import mbuild as mb
from mbuild.lib.atoms import H
from mbuild.lib.recipes import Monolayer
from mbuild.lib.surfaces import Betacristobalite
from mbuild.tests.base_test import BaseTest


class TestMonolayer(BaseTest):
    def test_monolayer(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)

        chain = mb.recipes.Polymer(monomers=[ch2])
        chain.build(n=10, add_hydrogens=False)
        monolayer = Monolayer(
            surface=Betacristobalite(),
            chains=chain,
            backfill=H(),
            pattern=pattern,
        )

        assert monolayer.n_particles == 2000 + n * m * 29
        assert monolayer.n_bonds == 2500 + n * m * 29

    def test_pattern_kwargs(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)

        chain = mb.recipes.Polymer(monomers=[ch2])
        chain.build(n=10, add_hydrogens=False)

        monolayer = Monolayer(
            surface=Betacristobalite(),
            chains=H(),
            guest_port_name="up",
            backfill=chain,
            backfill_port_name="down",
            pattern=pattern,
        )

        chains = 100 - (n * m)

        assert monolayer.n_particles == 2000 + chains * 29
        assert monolayer.n_bonds == 2500 + chains * 29

    def test_periodic_pattern(self, ch2):
        # Make Periodic without Hydrogen Conflict
        for axis in ["x", "y", "z"]:
            chain = mb.recipes.Polymer(monomers=[ch2])
            chain.build(n=10, add_hydrogens=False)
            chain.create_periodic_bond(axis=axis)
            assert not chain.all_ports()

        bonded_atoms = [x.name for x in list(chain["monomer[0]"][0].direct_bonds())]
        assert bonded_atoms.count("H") == 2
        assert bonded_atoms.count("C") == 2

        # Make Periodic with Hydrogen Conflict
        chain2 = mb.recipes.Polymer(monomers=[ch2])
        chain2.build(n=10, add_hydrogens=True)
        with pytest.raises(ValueError):
            chain2.create_periodic_bond(axis="y")

        # Make Periodic with End-Group Conflict
        chain3 = mb.recipes.Polymer(
            monomers=[ch2], end_groups=[mb.clone(ch2), mb.clone(ch2)]
        )
        chain3.build(n=10)
        with pytest.raises(ValueError):
            chain3.create_periodic_bond(axis="z")

        # Make Periodic with Unsupported Axis
        chain2 = mb.recipes.Polymer(monomers=[ch2])
        chain2.build(n=10, add_hydrogens=True)
        with pytest.raises(ValueError):
            chain2.create_periodic_bond(axis="a")

    def test_mixed_monolayer(self, ch2):
        n = 8
        m = 8
        pattern = mb.Grid2DPattern(n, m)
        fractions = [0.75, 0.25]

        chain_a = mb.recipes.Polymer(monomers=[ch2])
        chain_a.build(n=5, add_hydrogens=False)

        chain_b = mb.recipes.Polymer(monomers=[ch2])
        chain_b.build(n=15, add_hydrogens=False)

        monolayer = Monolayer(
            surface=Betacristobalite(),
            chains=[chain_a, chain_b],
            fractions=fractions,
            backfill=H(),
            pattern=pattern,
        )

        n_a = round(n * m * 0.75)
        n_b = round(n * m * 0.25)
        assert monolayer.n_particles == 2000 + n_a * 14 + n_b * 44
        assert monolayer.n_bonds == 2500 + n_a * 14 + n_b * 44
