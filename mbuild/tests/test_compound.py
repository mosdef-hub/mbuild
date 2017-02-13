import os

import numpy as np
import pytest

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.utils.io import get_fn, has_intermol
from mbuild.tests.base_test import BaseTest


class TestCompound(BaseTest):

    def test_load_and_create(self):
        mb.load(get_fn('methyl.pdb'))

    def test_update_from_file(self, ch3):
        ch3.update_coordinates(get_fn("methyl.pdb"))

    def test_save(self):
        methyl = mb.load(get_fn('methyl.pdb'))
        extensions = ['.xyz', '.pdb', '.mol2']
        for ext in extensions:
            outfile = 'methyl_out' + ext
            methyl.save(filename=outfile)
            assert os.path.exists(outfile)

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
        assert compound.n_particles == 8 + 2 * 3
        assert compound.n_bonds == 7 + 2 * 2

    def test_init_with_bad_name(self):
        with pytest.raises(ValueError):
            mb.Compound(name=1)

    def test_add_wrong_input(self, ethane):
        with pytest.raises(ValueError):
            ethane.add('water')

    def test_add_existing_parent(self, ethane, h2o):
        water_in_water = mb.clone(h2o)
        h2o.add(water_in_water)
        with pytest.raises(MBuildError):
            ethane.add(water_in_water)

    def test_add_label_exists(self, ethane, h2o):
        ethane.add(h2o, label='water')
        with pytest.raises(MBuildError):
            ethane.add(mb.clone(h2o), label='water')

    def test_set_pos(self, ethane):
        with pytest.raises(MBuildError):
            ethane.pos = [0, 0, 0]

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

    def test_generate_bonds(self, ch3):
        ch3.generate_bonds('H', 'H', dmin=0.01, dmax=2.0)
        assert ch3.n_bonds == 3 + 3

    def test_remove(self, ethane):
        hydrogens = ethane.particles_by_name('H')
        ethane.remove(hydrogens)

        assert ethane.n_particles == 2
        assert ethane.n_bonds == 1
        for part in ethane.children:
            assert part.n_bonds == 0

    def test_remove_no_bond_graph(self):
        compound = mb.Compound()
        particle = mb.Compound(name='C', pos=[0, 0, 0])
        compound.add(particle, 'test-particle')
        compound.remove(particle)
        assert particle not in compound.particles()

    def test_remove_bond(self, ch3):
        ch_bond = list(ch3.bonds())[0]
        ch3.remove_bond(ch_bond)
        assert ch3.n_bonds == 2

        with pytest.warns(UserWarning):
            ch3.remove_bond(ch_bond)

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

    def test_name(self):
        with pytest.raises(ValueError):
            mb.Compound(name=1)

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

    def test_reload(self):
        from mbuild.examples.pmpc.brush import Brush
        from numpy import pi
        # Create a compound and write it to file.
        brush1 = Brush()
        brush1.save("brush1.pdb")

        # Create another compound, rotate it and write it to file.
        brush2 = Brush()
        mb.rotate_around_z(brush2, pi/2)
        brush2.save("brush2.pdb")

        # Load brush2.pdb into brush1, modifying the atom positions of brush1.
        brush1.update_coordinates("brush2.pdb")
        brush1.save("modified_brush1.pdb")

        assert brush1['pmpc'].n_particles == 164
        assert brush1['pmpc'].n_bonds == 163
        assert len(brush1['pmpc']['monomer']) == 4
        assert brush1['pmpc']['monomer'][0].n_particles == 41
        assert brush1['pmpc']['monomer'][0].n_bonds == 40

    def test_to_trajectory(self, ethane, ch3):
        traj = ethane.to_trajectory()
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 1

        traj = ethane.to_trajectory(residue_types=ch3.__class__)
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

    @pytest.mark.skipif(not has_intermol, reason="InterMol is not installed")
    def test_intermol_conversion1(self, ethane, h2o):
        compound = mb.Compound([ethane, h2o])

        intermol_system = compound.to_intermol()
        assert len(intermol_system.molecule_types) == 1
        assert 'Compound' in intermol_system.molecule_types
        assert len(intermol_system.molecule_types['Compound'].bonds) == 9

        assert len(intermol_system.molecule_types['Compound'].molecules) == 1
        molecules = list(intermol_system.molecule_types['Compound'].molecules)
        assert len(molecules[0].atoms) == 11

    @pytest.mark.skipif(not has_intermol, reason="InterMol is not installed")
    def test_intermol_conversion2(self, ethane, h2o):
        # 2 distinct Ethane objects.
        compound = mb.Compound([ethane, mb.clone(ethane), h2o])

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

        # from_parmed tests
        compound2 = mb.Compound()
        compound2.from_parmed(structure)

        assert compound2.n_particles == 11
        assert len([at for at in compound2.particles() if at.name == 'C']) == 2
        assert len([at for at in compound2.particles() if at.name == 'H']) == 8
        assert len([at for at in compound2.particles() if at.name == 'O']) == 1

        assert compound2.n_bonds == 9

    def test_parmed_element_guess(self):
        compound = mb.Particle(name='foobar')
        with pytest.warns(UserWarning):
            _ = compound.to_parmed()

        compound = mb.Particle(name='XXXXXX')
        with pytest.warns(UserWarning):
            _ = compound.to_parmed()

    def test_min_periodic_dist(self, ethane):
        compound = mb.Compound(ethane)
        C_pos = np.array([atom.pos for atom in list(compound.particles_by_name('C'))])
        assert round(compound.min_periodic_distance(C_pos[0], C_pos[1]), 2) == 0.14
        compound.periodicity = np.array([0.2, 0.2, 0.2])
        assert round(compound.min_periodic_distance(C_pos[0], C_pos[1]), 2) == 0.06

    def test_bond_graph(self, ch3):
        compound = mb.Compound()
        compound.add(ch3)
        assert compound.n_bonds == 3
        assert all(compound.bond_graph.has_node(particle)
                   for particle in ch3.particles())

        ch3_nobonds = mb.clone(ch3)
        for bond in ch3_nobonds.bonds():
            ch3_nobonds.remove_bond(bond)
        compound.add(ch3_nobonds)
        assert compound.n_bonds == 3
        assert not any(compound.bond_graph.has_node(particle)
                       for particle in ch3_nobonds.particles())

        carbons = list(compound.particles_by_name('C'))
        compound.add_bond((carbons[0], carbons[1]))
        assert compound.n_bonds == 4
        assert all(compound.bond_graph.has_node(particle)
                   for particle in carbons)
        assert any(compound.bond_graph.has_node(particle)
                   for particle in ch3_nobonds.particles())

        compound.remove_bond((carbons[0], carbons[1]))
        assert not any(compound.bond_graph.has_node(particle)
                       for particle in ch3_nobonds.particles())
