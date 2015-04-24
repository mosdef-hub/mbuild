import sys

import numpy as np
import pytest

import mbuild as mb
from mbuild.components.small_groups.ch3 import CH3
from mbuild.examples.ethane.ethane import Ethane
from mbuild.utils.io import get_fn
from mbuild.tests.base_test import BaseTest


class TestCompound(BaseTest):

    def test_load_and_create(self):
        mb.load(get_fn('methyl.pdb'))

    def test_update_from_file(self):
        methyl = CH3()
        methyl.update_coordinates(get_fn("methyl.pdb"))

    def test_save(self):
        methyl = mb.load(get_fn('methyl.pdb'))
        methyl.save(filename='methyl_out.pdb')

    def test_batch_add(self, ethane, h2o):
        compound = mb.Compound()
        compound.add([ethane, h2o])
        assert compound.n_atoms == 8 + 3
        assert compound.n_bonds == 7 + 2

    @pytest.skipif(sys.version_info < (3, 0), reason='InterMol requires 2.7')
    def test_intermol_conversion1(self, ethane, h2o):
        compound = mb.Compound()
        compound.add([ethane, h2o])
        intermol_system = compound._to_intermol()
        assert len(intermol_system.molecule_types) == 1
        assert 'Compound' in intermol_system.molecule_types
        assert len(intermol_system.molecule_types['Compound'].bond_forces) == 9

        assert len(intermol_system.molecule_types['Compound'].molecules) == 1
        molecules = list(intermol_system.molecule_types['Compound'].molecules)
        assert len(molecules[0].atoms) == 11

    @pytest.skipif(sys.version_info < (3, 0), reason='InterMol requires 2.7')
    def test_intermol_conversion2(self, ethane, h2o):
        compound = mb.Compound()
        compound.add([ethane, Ethane(), h2o])  # 2 distinct Ethane objects
        molecule_types = [type(ethane), type(h2o)]
        intermol_system = compound._to_intermol(molecule_types=molecule_types)
        assert len(intermol_system.molecule_types) == 2
        assert 'Ethane' in intermol_system.molecule_types
        assert 'H2O' in intermol_system.molecule_types
        assert len(intermol_system.molecule_types['Ethane'].bond_forces) == 7
        assert len(intermol_system.molecule_types['H2O'].bond_forces) == 2

        assert len(intermol_system.molecule_types['Ethane'].molecules) == 2
        ethanes = list(intermol_system.molecule_types['Ethane'].molecules)
        assert len(ethanes[0].atoms) == len(ethanes[1].atoms) == 8

        assert len(intermol_system.molecule_types['H2O'].molecules) == 1
        h2os = list(intermol_system.molecule_types['H2O'].molecules)
        assert len(h2os[0].atoms) == 3

    def test_atom_list_by_kind(self, ethane):
        non_ports = ethane.atom_list_by_name()
        assert sum([1 for x in non_ports if x.name != 'G']) == 8

        with_G = ethane.atom_list_by_name(exclude_ports=False)
        assert len(with_G) == 24

        only_H = ethane.atom_list_by_name('H')
        assert sum([1 for _ in only_H]) == 6

        only_G = ethane.atom_list_by_name('G', exclude_ports=False)
        assert sum([1 for _ in only_G]) == 16

    def test_bond_list_by_kind(self, ethane):
        C_H_bonds = ethane.bond_list_by_kind(kind='C-H')
        assert sum([1 for x in C_H_bonds if x.kind == 'C-H']) == 6

        all_bonds = ethane.bond_list_by_kind()
        assert len(all_bonds) == 7

    # def test_atoms_in_range(self, ethane):
    #     group = ethane.atoms_in_range(ethane.atoms[0].pos, 0.141)
    #     assert sum([1 for x in group if x.name == 'H']) == 3
    #     assert sum([1 for x in group if x.name == 'C']) == 2
    #
    #     group = ethane.atoms_in_range(ethane.atoms[0].pos, 0.141, max_items=4)
    #     assert sum([1 for x in group if x.name == 'H']) == 3
    #     assert sum([1 for x in group if x.name == 'C']) == 1

    def test_remove(self, ethane):
        hydrogens = ethane.atom_list_by_name("H")
        ethane.remove(hydrogens)

        assert ethane.n_atoms == 2
        assert ethane.n_bonds == 1
        for atom in ethane.atoms:
            assert atom.n_bonds == 1

    def test_center(self, methane):
        assert np.array_equal(methane.center, np.array([0, 0, 0]))
        port = mb.Port()
        assert np.allclose(port.center, np.array([0.0, 0.0, 2.5e-3]))

    def test_visualize(self, ethane):
        ethane.visualize()

    def test_visualize_ports(self, ethane):
        ethane.visualize(show_ports=True)