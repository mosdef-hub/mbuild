from collections import Counter
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.validation import assert_port_exists

class TestPolymer(BaseTest):

    def test_polymer(self, ch2):
        n = 6
        c6 = mb.Polymer(ch2, n=n)
        assert c6.n_particles == n * 3
        assert c6.n_bonds == n * 2 + (n - 1)

    def test_block_copolymer(self, ch2, ester):
        n = 2
        sequence = 'ABBA'
        abba = mb.Polymer([ch2, ester], sequence=sequence, n=n)

        assert abba.n_particles == n * 3 * len(sequence)
        assert len(abba.children) == len(sequence) * n
        assert abba.children[0].name == 'CH2'
        assert abba.children[1].name == 'Ester'
        assert abba.children[2].name == 'Ester'
        assert abba.children[3].name == 'CH2'
        assert abba.children[4].name == 'CH2'
        assert abba.children[5].name == 'Ester'
        assert abba.children[6].name == 'Ester'
        assert abba.children[7].name == 'CH2'
        n_elements = Counter(p.name for p in abba.particles())
        assert n_elements['C'] == n * len(sequence)
        assert n_elements['H'] == n * len(sequence)
        assert n_elements['O'] == n * len(sequence)
        assert abba.n_bonds == n * 2 * len(sequence) + (n * len(sequence) - 1)

    def test_polymer_caps(self, ch2, ch3, ester):
        n = 6
        c6 = mb.Polymer(ch2, n=n,caps=(ch3,ester),cap_ports=('up','down'))
        assert c6.n_particles == n * 3 + 7
        assert c6.n_bonds == n * 2 + (n - 1) + 7
        assert c6.children[6].name == 'CH3'
        assert c6.children[7].name == 'Ester'
        assert len(c6.all_ports()) == 1
        assert assert_port_exists('up', c6.children[7])
        
        with pytest.raises(ValueError):
            assert_port_exists('up', c6.children[6])

    def test_block_copolymer_cap(self, ch2, ester, ch3):
        n = 2
        sequence = 'ABBA'
        abba = mb.Polymer([ch2, ester], sequence=sequence, n=n, caps = (None,ch3), cap_ports=(None,'up'))

        assert abba.n_particles == n * 3 * len(sequence) + 4
        assert len(abba.children) == len(sequence) * n + 1
        assert abba.children[0].name == 'CH2'
        assert abba.children[1].name == 'Ester'
        assert abba.children[2].name == 'Ester'
        assert abba.children[3].name == 'CH2'
        assert abba.children[4].name == 'CH2'
        assert abba.children[5].name == 'Ester'
        assert abba.children[6].name == 'Ester'
        assert abba.children[7].name == 'CH2'
        n_elements = Counter(p.name for p in abba.particles())
        assert n_elements['C'] == n * len(sequence) + 1
        assert n_elements['H'] == n * len(sequence) + 3
        assert n_elements['O'] == n * len(sequence)
        assert abba.n_bonds == n * 2 * len(sequence) + (n * len(sequence) - 1) + 4
        assert len(abba.all_ports()) == 1
        assert assert_port_exists('down',abba)

    def test_cap_values(self, ch2, ch3):
        n=2
        with pytest.raises(ValueError):
            mb.Polymer(ch2, n=n, caps = (None,ch3), cap_ports=('up',None))