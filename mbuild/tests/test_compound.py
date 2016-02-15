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
        assert compound.n_particles == 8 + 3
        assert compound.n_bonds == 7 + 2

    def test_init_with_subcompounds1(self, ethane):
        compound = mb.Compound(ethane)
        assert compound.n_particles == 8
        assert compound.n_bonds == 7

    def test_init_with_subcompounds2(self, ethane, h2o):
        compound = mb.Compound([ethane, h2o])
        assert compound.n_particles == 8 + 3
        assert compound.n_bonds == 7 + 2

    def test_init_with_subcompounds3(self, ethane, h2o):
        compound = mb.Compound([ethane, [h2o, mb.clone(h2o)]])
        assert compound.n_particles == 8 + 2*3
        assert compound.n_bonds == 7 + 2 * 2

    def test_xyz(self, ethane):
        xyz = ethane.xyz
        assert xyz.shape == (8, 3)

        xyz = ethane.xyz_with_ports
        assert xyz.shape == (24, 3)

    def test_particles_by_name(self, ethane):
        assert sum(1 for _ in ethane.particles()) == 8

        only_H = ethane.particles_by_name('H')
        assert sum(1 for _ in only_H) == 6

        only_C = ethane.particles_by_name('C')
        assert sum(1 for _ in only_C) == 2

    def test_particles_in_range(self, ethane):
        group = ethane.particles_in_range(ethane[0], 0.141)
        assert sum([1 for x in group if x.name == 'H']) == 3
        assert sum([1 for x in group if x.name == 'C']) == 2

        group = ethane.particles_in_range(ethane[0], 0.141, max_particles=4)
        assert sum([1 for x in group if x.name == 'H']) == 3
        assert sum([1 for x in group if x.name == 'C']) == 1

    def test_add_bonds(self, ch3):
        ch3.generate_bonds('H', 'H', dmin=0.01, dmax=2.0)
        assert ch3.n_bonds == 3 + 3

    def test_remove(self, ethane):
        hydrogens = ethane.particles_by_name('H')
        ethane.remove(hydrogens)

        assert ethane.n_particles == 2
        assert ethane.n_bonds == 1
        for part in ethane.children:
            assert part.n_bonds == 0

    def test_center(self, methane):
        assert np.array_equal(methane.center, np.array([0, 0, 0]))
        port = mb.Port()
        assert np.allclose(port.center, np.array([0.0, 0.0, 2.5e-3]))

    def test_single_particle(self):
        part = mb.Particle(name='A')
        assert part.n_particles == 1
        assert len(list(part.particles())) == 1
        assert part.xyz.shape == (1, 3)
        assert part.root == part
        assert len(list(part.ancestors())) == 0
        assert next(part.particles_by_name('A')) == part

    def test_particle_in_particle(self):
        part = mb.Particle(name='A')
        parent = mb.Compound(part)

        assert part.n_particles == 1
        assert len(list(part.particles())) == 1
        assert part.xyz.shape == (1, 3)
        assert part.root == parent
        assert len(list(part.ancestors())) == 1
        assert next(part.particles_by_name('A')) == part

        assert parent.n_particles == 1
        assert len(list(parent.particles())) == 1
        assert parent.xyz.shape == (1, 3)
        assert parent.root == parent
        assert len(list(parent.ancestors())) == 0
        assert next(parent.particles_by_name('A')) == part

    @pytest.mark.skipif(bool(os.getenv("CI")), reason="Running on CI")
    def test_visualize(self, ethane):
        ethane.visualize()

    @pytest.mark.skipif(bool(os.getenv("CI")), reason="Running on CI")
    def test_visualize_ports(self, ethane):
        ethane.visualize(show_ports=True)

    # Conversions
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

        methyl = next(iter(ethane.children))
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

    @pytest.mark.skipif(bool(os.getenv("CI")), reason="Running on CI")
    def test_intermol_conversion1(self, ethane, h2o):
        compound = mb.Compound([ethane, h2o])

        intermol_system = compound.to_intermol()
        assert len(intermol_system.molecule_types) == 1
        assert 'Compound' in intermol_system.molecule_types
        assert len(intermol_system.molecule_types['Compound'].bonds) == 9

        assert len(intermol_system.molecule_types['Compound'].molecules) == 1
        molecules = list(intermol_system.molecule_types['Compound'].molecules)
        assert len(molecules[0].atoms) == 11

    @pytest.mark.skipif(bool(os.getenv("CI")), reason="Running on CI")
    def test_intermol_conversion2(self, ethane, h2o):
        compound = mb.Compound([ethane, mb.clone(ethane), h2o]) # 2 distinct Ethane objects

        molecule_types = [type(ethane), type(h2o)]
        intermol_system = compound.to_intermol(molecule_types=molecule_types)
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

    def test_parmed_conversion(self, ethane, h2o):
        compound = mb.Compound([ethane, h2o])

        structure = compound.to_parmed()
        assert structure.title == 'Compound'

        structure = compound.to_parmed(title='eth_h2o')
        assert structure.title == 'eth_h2o'

        assert len(structure.atoms) == 11
        assert len([at for at in structure.atoms if at.element == 6]) == 2
        assert len([at for at in structure.atoms if at.element == 1]) == 8
        assert len([at for at in structure.atoms if at.element == 8]) == 1

        assert len(structure.bonds) == 9

    def test_parmed_element_guess(self):
        compound = mb.Particle(name='foobar')
        with pytest.warns(UserWarning):
            structure = compound.to_parmed()

        compound = mb.Particle(name='XXXXXX')
        with pytest.warns(UserWarning):
            structure = compound.to_parmed()

    def test_min_periodic_dist(self,ethane):
        compound = mb.Compound(ethane)
        C_pos = np.array([atom.pos for atom in list(compound.particles()) if atom.name == 'C'])
        assert round(compound.min_periodic_distance(C_pos[0],C_pos[1]),2) == 0.14
        compound.periodicity = np.array([0.2,0.2,0.2])
        assert round(compound.min_periodic_distance(C_pos[0],C_pos[1]),2) == 0.06
