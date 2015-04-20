import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestPolymer(BaseTest):

    def test_polymer(self, ch2):
        n = 6
        c6 = mb.Polymer(ch2, n=n)
        assert c6.n_atoms == n * 3
        assert c6.n_bonds == n * 2 + (n - 1)
