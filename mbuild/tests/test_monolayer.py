import mbuild as mb
from mbuild.tests.base_test import BaseTest

from mbuild.components.surfaces.betacristobalite import Betacristobalite
from mbuild.components.atoms.H import H


class TestMonolayer(BaseTest):

    def test_monolayer(self, ch2):
        n = 8
        m = 8
        mask = mb.grid_mask_2d(n, m)

        chain = mb.Polymer(ch2, n=10)
        print(chain.n_bonds)
        monolayer = mb.Monolayer(surface=Betacristobalite(), chain=chain,
                                 backfill=H(), mask=mask)

        assert monolayer.n_atoms == 1800 + n * m * (10*3) + (100 - n*m)
        assert monolayer.n_bonds == 2300 + n * m * (10*2 + 9 + 1) + (100 - n*m)
