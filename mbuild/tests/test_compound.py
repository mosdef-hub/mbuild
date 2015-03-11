import pytest

import mbuild as mb
from mbuild.components.small_groups.ch3 import Ch3
from mbuild.testing.tools import get_fn
from mbuild.tests.base_test import BaseTest


class TestCompound(BaseTest):

    def test_load_and_create(self):
        mb.load(get_fn('methyl.pdb'))

    def test_update_from_file(self):
        methyl = Ch3()
        methyl.update_coordinates(get_fn("methyl.pdb"))

    def test_save(self):
        methyl = mb.load(get_fn('methyl.pdb'))
        methyl.save(filename='methyl_out.pdb')

    @pytest.fixture
    def ethane(self):
        from mbuild.examples.ethane.ethane import Ethane
        return Ethane()

    def test_atom_list_by_kind(self, ethane):
        non_ports = ethane.atom_list_by_kind(excludeG=True)
        assert sum([1 for x in non_ports if x.kind != 'G']) == 8

        with_G = ethane.atom_list_by_kind()
        assert len(with_G) == 24

        only_H = ethane.atom_list_by_kind('H')
        assert sum([1 for x in only_H if x.kind == 'H']) == 6

        only_G = ethane.atom_list_by_kind('G')
        assert sum([1 for x in only_G if x.kind == 'G']) == 16

    def test_bond_list_by_kind(self, ethane):
        C_H_bonds = ethane.bond_list_by_kind(kind='C-H')
        assert sum([1 for x in C_H_bonds if x.kind == 'C-H']) == 6

        all_bonds = ethane.bond_list_by_kind()
        assert len(all_bonds) == 7

    def test_atoms_in_range(self, ethane):
        group = ethane.atoms_in_range(ethane.atoms[0].pos, 0.141)
        assert sum([1 for x in group if x.kind == 'H']) == 3
        assert sum([1 for x in group if x.kind == 'C']) == 2

        group = ethane.atoms_in_range(ethane.atoms[0].pos, 0.141, max_items=4)
        assert sum([1 for x in group if x.kind == 'H']) == 3
        assert sum([1 for x in group if x.kind == 'C']) == 1

    def test_remove(self, ethane):
        hydrogens = ethane.atom_list_by_kind("H")
        ethane.remove(hydrogens)

        assert ethane.n_atoms == 2
        assert ethane.n_bonds == 1
        for atom in ethane.atoms:
            assert atom.n_bonds == 1
