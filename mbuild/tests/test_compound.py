import os

import numpy as np
import parmed as pmd
import pytest

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.utils.io import get_fn, has_intermol, has_openbabel
from mbuild.tests.base_test import BaseTest

class TestCompound(BaseTest):

    def test_load_and_create(self):
        mb.load(get_fn('methyl.pdb'))

    def test_update_from_file(self, ch3):
        ch3.update_coordinates(get_fn("methyl.pdb"))

    def test_save_simple(self, ch3):
        extensions = ['.xyz', '.pdb', '.mol2']
        for ext in extensions:
            outfile = 'methyl_out' + ext
            ch3.save(filename=outfile)
            assert os.path.exists(outfile)

    def test_save_box(self, ch3):
        extensions = ['.mol2', '.pdb', '.hoomdxml', '.gro']
        box_attributes = ['mins', 'maxs', 'lengths']
        custom_box = mb.Box([.8, .8, .8])
        for ext in extensions:
            outfile_padded = 'padded_methyl' + ext
            outfile_custom = 'custom_methyl' + ext
            ch3.save(filename=outfile_padded, box=None, overwrite=True)
            ch3.save(filename=outfile_custom, box=custom_box, overwrite=True)
            padded_ch3 = mb.load(outfile_padded)
            custom_ch3 = mb.load(outfile_custom)
            for attr in box_attributes:
                pad_attr = getattr(padded_ch3.boundingbox, attr)
                custom_attr = getattr(custom_ch3.boundingbox, attr)
                assert np.array_equal(pad_attr, custom_attr)

    def test_save_overwrite(self, ch3):
        extensions = ['.gsd', '.hoomdxml', '.lammps', '.lmp', '.top', '.gro']
        for ext in extensions:
            outfile = 'lyhtem' + ext
            ch3.save(filename=outfile)
            ch3.save(filename=outfile, overwrite=True)
            with pytest.raises(IOError):
                ch3.save(filename=outfile, overwrite=False)

    def test_save_forcefield(self, methane):
        exts = ['.gsd', '.hoomdxml', '.lammps', '.lmp', '.top', '.gro',
                '.mol2', '.pdb', '.xyz']
        for ext in exts:
            methane.save('lythem' + ext,
                         forcefield_name='oplsaa',
                         overwrite=True)

    def test_save_resnames(self, ch3, h2o):
        system = mb.Compound([ch3, h2o])
        system.save('resnames.gro', residues=['CH3', 'H2O'])
        struct = pmd.load_file('resnames.gro')

        assert struct.residues[0].name == 'CH3'
        assert struct.residues[1].name == 'H2O'

    def test_save_resnames_single(self, c3, n4):
        system = mb.Compound([c3, n4])
        system.save('resnames_single.gro', residues=['C3', 'N4'])
        struct = pmd.load_file('resnames_single.gro')
        assert struct.residues[0].number ==  1
        assert struct.residues[1].number ==  2

    def test_save_references(self, methane):
        methane.save('methyl.mol2', forcefield_name='oplsaa',
                     references_file='methane.bib')
        assert os.path.isfile('methane.bib')

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

    def test_remove_from_box(self, ethane):
        n_ethanes = 5
        box = mb.fill_box(ethane, n_ethanes, [3, 3, 3])
        box.remove(box.children[3])

        n_ethanes -= 1
        assert box.n_particles == n_ethanes * ethane.n_particles
        assert len(box.children) == n_ethanes
        assert box.n_bonds == n_ethanes * ethane.n_bonds
        assert len([meth.referenced_ports()
                    for eth in box.children
                    for meth in eth.children]) == 2 * n_ethanes

    def test_remove(self, ethane):
        hydrogens = ethane.particles_by_name('H')
        ethane.remove(hydrogens)

        assert ethane.n_particles == 2
        assert ethane.n_bonds == 1
        for part in ethane.children:
            assert part.n_bonds == 0

        carbons = ethane.particles_by_name('C')
        ethane.remove(carbons)
        assert ethane.n_particles == 0
        assert ethane.n_bonds == 0
        assert len(ethane.children) == 2
        assert len(ethane.children[0].children) == 1  # Still contains a port

    def test_remove_many(self, ethane):
        ethane.remove([ethane.children[0], ethane.children[1]])

        assert ethane.n_particles == 1
        assert ethane._n_particles() == 0
        assert ethane.n_bonds == 0
        for part in ethane.children:
            assert isinstance(part, mb.Port)

    def test_remove_subcompound(self, ethane):
        methyl = ethane.children[0]
        ethane.remove(methyl)

        assert ethane.n_particles == 4
        assert ethane.n_bonds == 3
        assert len(ethane.children) == 1
        assert len(ethane.children[0].children) == 5  # Still contains a port

        methyl = ethane.children[0]
        ethane.remove(methyl)

        assert ethane.n_particles == 1
        assert ethane._n_particles() == 0
        assert ethane.n_bonds == 0
        assert len(ethane.children) == 0

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
        for orientation in np.identity(3):
            separation = 0.2
            port = mb.Port(anchor=methane[0], orientation=orientation)
            assert np.allclose(port.center, np.array([0.0, 0.0, 0.0]), atol=1e-15)
            port = mb.Port(anchor=methane[0], orientation=orientation,
                           separation=separation)
            assert np.allclose(port.center, separation*orientation, atol=1e-15)
        np.random.seed(0)
        for orientation in np.random.rand(5, 3):
            port = mb.Port(anchor=methane[0], orientation=orientation)
            assert np.allclose(port.center, np.array([0.0, 0.0, 0.0]), atol=1e-15)
            port = mb.Port(anchor=methane[0], orientation=orientation,
                           separation=separation)
            assert np.allclose(port.center,
                               separation*orientation/np.linalg.norm(orientation),
                               atol=1e-15)

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
        mb.rotate(brush2, pi/2, [0, 0, 1])
        brush2.save("brush2.pdb")

        # Load brush2.pdb into brush1, modifying the atom positions of brush1.
        brush1.update_coordinates("brush2.pdb")
        brush1.save("modified_brush1.pdb")

        assert brush1['pmpc'].n_particles == 164
        assert brush1['pmpc'].n_bonds == 163
        assert len(brush1['pmpc']['monomer']) == 4
        assert brush1['pmpc']['monomer'][0].n_particles == 41
        assert brush1['pmpc']['monomer'][0].n_bonds == 40

    def test_to_trajectory(self, ethane, c3, n4):
        traj = ethane.to_trajectory()
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 1

        traj = ethane.to_trajectory(residues='CH3')
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 2
        assert 'CH3' in [res.name for res in traj.top.residues]
        assert all(res.n_atoms == 4 for res in traj.top.residues)

        traj = ethane.to_trajectory(chains='CH3')
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 2
        assert traj.n_residues == 2
        assert all(chain.n_atoms == 4 for chain in traj.top.chains)
        assert all(chain.n_residues == 1 for chain in traj.top.chains)

        system = mb.Compound([c3, n4])
        traj = system.to_trajectory(residues=['C', 'N'])
        assert traj.n_atoms == 2
        assert traj.top.n_bonds == 0
        assert traj.n_chains == 1
        assert traj.n_residues == 2

        traj = system.to_trajectory(chains=['C', 'N'])
        assert traj.n_atoms == 2
        assert traj.top.n_bonds == 0
        assert traj.n_chains == 2
        assert traj.n_residues == 2

        methyl = next(iter(ethane.children))
        traj = methyl.to_trajectory()
        assert traj.n_atoms == 4
        assert traj.top.n_bonds == 3
        assert traj.n_chains == 1
        assert traj.n_residues == 1

    def test_resnames_mdtraj(self, h2o, ethane):
        system = mb.Compound([h2o, mb.clone(h2o), ethane])
        traj = system.to_trajectory(residues=['Ethane', 'H2O'])
        residues = list(traj.top.residues)
        assert traj.n_residues == 3
        assert residues[0].name == 'H2O'
        assert residues[1].name == 'H2O'
        assert residues[2].name == 'Ethane'

        traj = system.to_trajectory(residues='Ethane')
        residues = list(traj.top.residues)
        assert traj.n_residues == 2
        assert residues[0].name == 'RES'
        assert residues[1].name == 'Ethane'

        traj = system.to_trajectory(residues=['Ethane'])
        residues = list(traj.top.residues)
        assert traj.n_residues == 2
        assert residues[0].name == 'RES'
        assert residues[1].name == 'Ethane'

        traj = system.to_trajectory()
        residues = list(traj.top.residues)
        assert traj.n_residues == 1
        assert residues[0].name == 'RES'

    def test_chainnames_mdtraj(self, h2o, ethane):
        system = mb.Compound([h2o, mb.clone(h2o), ethane])
        traj = system.to_trajectory(chains=['Ethane', 'H2O'])
        assert traj.n_chains == 3

        traj = system.to_trajectory(chains='Ethane')
        assert traj.n_chains == 2

        traj = system.to_trajectory(chains=['Ethane'])
        assert traj.n_chains == 2

        traj = system.to_trajectory()
        assert traj.n_chains == 1

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

        assert (sum(len(res.atoms) for res in structure.residues) ==
                len(structure.atoms))

        compound2 = mb.Compound()
        compound2.from_parmed(structure)

        assert compound2.n_particles == 11
        assert len([at for at in compound2.particles() if at.name == 'C']) == 2
        assert len([at for at in compound2.particles() if at.name == 'H']) == 8
        assert len([at for at in compound2.particles() if at.name == 'O']) == 1

        assert compound2.n_bonds == 9

    def test_resnames_parmed(self, h2o, ethane):
        system = mb.Compound([h2o, mb.clone(h2o), ethane])
        struct = system.to_parmed(residues=['Ethane', 'H2O'])
        assert len(struct.residues) == 3
        assert struct.residues[0].name == 'H2O'
        assert struct.residues[1].name == 'H2O'
        assert struct.residues[2].name == 'Ethane'
        assert sum(len(res.atoms) for res in struct.residues) == len(struct.atoms)

        struct = system.to_parmed(residues=['Ethane', 'H2O'])
        assert len(struct.residues) == 3
        assert struct.residues[0].name == 'H2O'
        assert struct.residues[1].name == 'H2O'
        assert struct.residues[2].name == 'Ethane'
        assert sum(len(res.atoms) for res in struct.residues) == len(struct.atoms)

        struct = system.to_parmed(residues='Ethane')
        assert len(struct.residues) == 2
        assert struct.residues[0].name == 'RES'
        assert struct.residues[1].name == 'Ethane'
        assert sum(len(res.atoms) for res in struct.residues) == len(struct.atoms)

        struct = system.to_parmed()
        assert len(struct.residues) == 1
        assert struct.residues[0].name == 'RES'
        assert sum(len(res.atoms) for res in struct.residues) == len(struct.atoms)

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

    def test_update_coords_update_ports(self, ch2):
        distances = np.round([ch2.min_periodic_distance(port.pos, ch2[0].pos)
                              for port in ch2.referenced_ports()], 5)
        orientations = np.round([port.pos - port.anchor.pos
                                 for port in ch2.referenced_ports()], 5)

        ch2_clone = mb.clone(ch2)
        ch2_clone[0].pos += [1, 1, 1]
        ch2_clone.save('ch2-shift.pdb')

        ch2.update_coordinates('ch2-shift.pdb')
        updated_distances = np.round([ch2.min_periodic_distance(port.pos, ch2[0].pos)
                                      for port in ch2.referenced_ports()], 5)
        updated_orientations = np.round([port.pos - port.anchor.pos
                                         for port in ch2.referenced_ports()], 5)

        assert np.array_equal(distances, updated_distances)
        assert np.array_equal(orientations, updated_orientations)

    def test_charge(self, ch2, ch3):
        compound = mb.Compound(charge=2.0)
        assert compound.charge == 2.0
        compound2 = mb.Compound()
        assert compound2.charge == 0.0

        ch2[0].charge = 0.5
        ch2[1].charge = -0.25
        ch3[0].charge = 1.0
        compound.add([ch2, ch3])
        assert compound.charge == 1.25
        assert ch2.charge == 0.25
        assert compound[0].charge == 0.5

        with pytest.raises(AttributeError):
            compound.charge = 2.0

    def test_charge_subcompounds(self, ch2, ch3):
        ch2[0].charge = 0.5
        ch2[1].charge = -0.25
        compound = mb.Compound(subcompounds=ch2)
        assert compound.charge == 0.25

        with pytest.raises(MBuildError):
            compound = mb.Compound(subcompounds=ch3, charge=1.0)

    def test_charge_neutrality_warn(self, benzene):
        benzene[0].charge = 0.25
        with pytest.warns(UserWarning):
            benzene.save('charge-test.mol2')

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization(self, octane):
        octane.energy_minimization()

    @pytest.mark.skipif(has_openbabel, reason="Open Babel package is installed")
    def test_energy_minimization_openbabel_warn(self, octane):
        with pytest.raises(MBuildError):
            octane.energy_minimization()

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization_ff(self, octane):
        for ff in ['UFF', 'GAFF', 'MMFF94', 'MMFF94s', 'Ghemical']:
            octane.energy_minimization(forcefield=ff)
        with pytest.raises(MBuildError):
            octane.energy_minimization(forcefield='fakeFF')

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization_algorithm(self, octane):
        for algorithm in ['cg', 'steep', 'md']:
            octane.energy_minimization(algorithm=algorithm)
        with pytest.raises(MBuildError):
            octane.energy_minimization(algorithm='fakeAlg')

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization_non_element(self, octane):
        for particle in octane.particles():
            particle.name = 'Q'
        with pytest.raises(MBuildError):
            octane.energy_minimization()

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization_ports(self, octane):
        distances = np.round([octane.min_periodic_distance(port.pos, port.anchor.pos)
                              for port in octane.all_ports()], 5)
        orientations = np.round([port.pos - port.anchor.pos
                                 for port in octane.all_ports()], 5)

        octane.energy_minimization()

        updated_distances = np.round([octane.min_periodic_distance(port.pos,
                                                                   port.anchor.pos)
                                      for port in octane.all_ports()], 5)
        updated_orientations = np.round([port.pos - port.anchor.pos
                                         for port in octane.all_ports()], 5)

        assert np.array_equal(distances, updated_distances)
        assert np.array_equal(orientations, updated_orientations)

    def test_clone_outside_containment(self, ch2, ch3):
        compound = mb.Compound()
        compound.add(ch2)
        mb.force_overlap(ch3, ch3['up'], ch2['up'])
        with pytest.raises(MBuildError):
            ch3_clone = mb.clone(ch3)

    def test_load_mol2_mdtraj(self):
        with pytest.raises(KeyError):
            mb.load(get_fn('benzene-nonelement.mol2'))
        mb.load(get_fn('benzene-nonelement.mol2'), use_parmed=True)

    def test_siliane_bond_number(self, silane):
        assert silane.n_bonds == 4
