import os
from collections import Counter

import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestPolymer(BaseTest):
    def test_polymer_from_smiles(self):
        chain = mb.recipes.Polymer()
        ethane = mb.load("CC", smiles=True)
        chain.add_monomer(ethane, indices=[2, -2], separation=0.15, replace=True)
        chain.build(n=5, add_hydrogens=True)
        assert len([p for p in chain.particles() if p.name == "C"]) == 10
        assert len([p for p in chain.particles() if p.name == "H"]) == 22
        assert len(chain.available_ports()) == 0
        assert len(chain.children) == 5

    def test_add_end_groups(self, ch2, ester):
        n = 6
        c6 = mb.recipes.Polymer(monomers=[ch2])
        acid = mb.load("C(=O)O", smiles=True)
        c6.add_end_groups(acid, index=3, separation=0.15)
        c6.build(n=n, add_hydrogens=False)
        assert len([p for p in c6.particles() if p.name == "O"]) == 4

    def test_pass_end_groups(self, ch2, ester):
        ester_2 = mb.clone(ester)
        c6 = mb.recipes.Polymer(monomers=[ch2], end_groups=[ester, ester_2])
        c6.build(n=6)
        assert c6.children[-1].name == "Ester"
        assert c6.children[-2].name == "Ester"

    def test_errors(self, ch2, ester):
        with pytest.raises(ValueError):  # Not enough end groups
            chain = mb.recipes.Polymer(monomers=[ch2], end_groups=[ester])

        with pytest.raises(ValueError):  # Bad sequence
            chain = mb.recipes.Polymer(monomers=[ch2])
            chain.build(n=5, sequence="AB")

        with pytest.raises(ValueError):  # Bad n value
            chain = mb.recipes.Polymer(monomers=[ch2])
            chain.build(n=0, sequence="A")

        with pytest.raises(ValueError):  # Bad end group label
            chain = mb.recipes.Polymer(monomers=[ch2])
            acid = mb.load("C(=O)O", smiles=True)
            chain.add_end_groups(
                acid, index=3, separation=0.15, duplicate=False, label="front"
            )

    def test_no_end_groups(self):
        chain = mb.recipes.Polymer()
        ethane = mb.load("CC", smiles=True)
        chain.add_monomer(ethane, indices=[2, -2], separation=0.15, replace=True)
        chain.build(n=5, add_hydrogens=False)
        assert len([p for p in chain.particles() if p.name == "H"]) == 20
        assert len(chain.available_ports()) == 2

    def test_replace_is_false(self):
        n = 6
        ch2 = mb.load(os.path.join(mb.__path__[0], "lib/moieties/ch2.pdb"))
        chain = mb.recipes.Polymer()
        chain.add_monomer(
            ch2,
            indices=[0, 0],
            orientation=[[0, 1, 0], [0, -1, 0]],
            separation=0.15,
            replace=False,
        )
        chain.build(n=n, add_hydrogens=False)
        assert chain.n_particles == n * 3
        assert chain.n_bonds == n * 2 + (n - 1)

    def test_polymer_from_moieties(self, ch2):
        n = 6
        c6 = mb.recipes.Polymer(monomers=[ch2])
        c6.build(n=n, add_hydrogens=False)
        assert c6.n_particles == n * 3
        assert c6.n_bonds == n * 2 + (n - 1)

    def test_block_copolymer(self, ch2, ester):
        n = 2
        sequence = "ABBA"
        abba = mb.recipes.Polymer(monomers=[ch2, ester])
        abba.build(n=n, sequence=sequence, add_hydrogens=False)

        assert abba.n_particles == n * 3 * len(sequence)
        assert len(abba.children) == len(sequence) * n
        assert abba.children[0].name == "CH2"
        assert abba.children[1].name == "Ester"
        assert abba.children[2].name == "Ester"
        assert abba.children[3].name == "CH2"
        assert abba.children[4].name == "CH2"
        assert abba.children[5].name == "Ester"
        assert abba.children[6].name == "Ester"
        assert abba.children[7].name == "CH2"
        n_elements = Counter(p.name for p in abba.particles())
        assert n_elements["C"] == n * len(sequence)
        assert n_elements["H"] == n * len(sequence)
        assert n_elements["O"] == n * len(sequence)
        assert abba.n_bonds == n * 2 * len(sequence) + (n * len(sequence) - 1)
