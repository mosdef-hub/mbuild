import json
import os
import numpy as np
import pytest
import mbuild as mb
from mbuild.utils.io import get_fn
from mbuild.tests.base_test import BaseTest

class TestCompound(BaseTest):

    def test_load_and_create(self):
        mb.load(get_fn('methyl.pdb'))

    def test_update_from_file(self, ch3):
        ch3.update_coordinates(get_fn("methyl.pdb"))

    def test_save(self):
        methyl = mb.load(get_fn('methyl.pdb'))
        methyl.save(filename='methyl_out.pdb')

    def test_batch_add(self, ethane, h2o):
        compound = mb.Compound()
        compound.add([ethane, h2o])
        assert compound.n_atoms == 8 + 3
        assert compound.n_bonds == 7 + 2

    def test_init_with_subcompounds1(self, ethane):
        compound = mb.Compound(ethane)
        assert compound.n_atoms == 8
        assert compound.n_bonds == 7

    def test_init_with_subcompounds2(self, ethane, h2o):
        compound = mb.Compound([ethane, h2o])
        assert compound.n_atoms == 8 + 3
        assert compound.n_bonds == 7 + 2

    def test_init_with_subcompounds3(self, ethane, h2o):
        compound = mb.Compound([ethane, [h2o, mb.clone(h2o)]])
        assert compound.n_atoms == 8 + 2*3
        assert compound.n_bonds == 7 + 2*2

    @pytest.mark.skipif(True, reason='Waiting for InterMol to stabilize')
    def test_intermol_conversion1(self, ethane, h2o):
        compound = mb.Compound()
        compound.add([ethane, h2o])
        intermol_system = compound._to_intermol()
        assert len(intermol_system.molecule_types) == 1
        assert 'Compound' in intermol_system.molecule_types
        assert len(intermol_system.molecule_types['Compound'].bonds) == 9

        assert len(intermol_system.molecule_types['Compound'].molecules) == 1
        molecules = list(intermol_system.molecule_types['Compound'].molecules)
        assert len(molecules[0].atoms) == 11

    @pytest.mark.skipif(True, reason='Waiting for InterMol to stabilize')
    def test_intermol_conversion2(self, ethane, h2o):
        compound = mb.Compound()
        compound.add([ethane, mb.clone(ethane), h2o])  # 2 distinct Ethane objects
        molecule_types = [type(ethane), type(h2o)]
        intermol_system = compound._to_intermol(molecule_types=molecule_types)
        assert len(intermol_system.molecule_types) == 2
        assert 'Ethane' in intermol_system.molecule_types
        assert 'H2O' in intermol_system.molecule_types
        assert len(intermol_system.molecule_types['Ethane'].bonds) == 7
        assert len(intermol_system.molecule_types['H2O'].bonds) == 2

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

    def test_atoms_in_range(self, ethane):
        group = ethane.atoms_in_range(ethane.atoms[0], 0.141)
        assert sum([1 for x in group if x.name == 'H']) == 3
        assert sum([1 for x in group if x.name == 'C']) == 2

        group = ethane.atoms_in_range(ethane.atoms[0], 0.141, max_atoms=4)
        assert sum([1 for x in group if x.name == 'H']) == 3
        assert sum([1 for x in group if x.name == 'C']) == 1

    def test_add_bonds(self, ch3):
        ch3.add_bonds('H', 'H', dmin=0.01, dmax=2.0)
        assert ch3.n_bonds == 3 + 3

    def test_remove(self, ethane):
        hydrogens = ethane.atom_list_by_name('H')
        ethane.remove(hydrogens)

        assert ethane.n_atoms == 2
        assert ethane.n_bonds == 1
        for atom in ethane.atoms:
            assert atom.n_bonds == 1

    def test_center(self, methane):
        assert np.array_equal(methane.center, np.array([0, 0, 0]))
        port = mb.Port()
        assert np.allclose(port.center, np.array([0.0, 0.0, 2.5e-3]))

    @pytest.mark.skipif(bool(os.getenv("CI")), reason="Running on CI")
    def test_visualize(self, ethane):
        ethane.visualize()

    @pytest.mark.skipif(bool(os.getenv("CI")), reason="Running on CI")
    def test_visualize_ports(self, ethane):
        ethane.visualize(show_ports=True)

    def test_to_trajectory(self, ethane, ch3):
        traj = ethane.to_trajectory()
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 1

        traj = ethane.to_trajectory(residue_types=ch3)
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 2
        assert 'CH3' in [res.name for res in traj.top.residues]
        assert all(res.n_atoms == 4 for res in traj.top.residues)

        traj = ethane.to_trajectory(chain_types=ch3)
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 2
        assert traj.n_residues == 2
        assert all(chain.n_atoms == 4 for chain in traj.top.chains)
        assert all(chain.n_residues == 1 for chain in traj.top.chains)

        methyl = next(iter(ethane.parts))
        traj = methyl.to_trajectory()
        assert traj.n_atoms == 4
        assert traj.top.n_bonds == 3
        assert traj.n_chains == 1
        assert traj.n_residues == 1

    def test_to_json(self, ethane):
        output = json.loads(ethane._to_json())
        assert len(output['atoms']) == 8
        assert len(output['bonds']) == 7

        output = json.loads(ethane._to_json(show_ports=True))
        assert len(output['atoms']) == 8+16
        assert len(output['bonds']) == 7
