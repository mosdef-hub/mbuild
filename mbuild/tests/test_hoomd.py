import numpy as np
import pytest
import xml.etree.ElementTree

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, has_foyer, has_hoomd, import_


@pytest.mark.skipif(not has_hoomd, reason="HOOMD is not installed")
class TestHoomd(BaseTest):
    def test_compound_to_snapshot(self, ethane):
        hoomd_snapshot = import_("mbuild.formats.hoomd_snapshot")
        snap, _ = hoomd_snapshot.to_hoomdsnapshot(ethane)

        assert snap.particles.N == 8
        assert snap.bonds.N == 7
        assert snap.angles.N == 0

    def test_particles_to_snapshot(self):
        hoomd_snapshot = import_("mbuild.formats.hoomd_snapshot")
        part = mb.Compound(name='Ar')
        box = mb.fill_box(part, n_compounds=10, box=mb.Box([5,5,5]))
        snap, _ = hoomd_snapshot.to_hoomdsnapshot(box)

        assert snap.particles.N == 10
        assert snap.bonds.N == 0
        assert snap.angles.N == 0


    def test_bad_input_to_snapshot(self):
        hoomd_snapshot = import_("mbuild.formats.hoomd_snapshot")
        with pytest.raises(ValueError):
            hoomd_snapshot.to_hoomdsnapshot('fake_object')

    def test_non_param_struc_to_snapshot(self, ethane):
        hoomd_snapshot = import_("mbuild.formats.hoomd_snapshot")
        structure = ethane.to_parmed()
        snap,_ = hoomd_snapshot.to_hoomdsnapshot(structure)

        assert snap.particles.N == 8
        assert snap.bonds.N == 7
        assert snap.angles.N == 0

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_param_structure_to_snapshot(self, ethane):
        hoomd_snapshot = import_("mbuild.formats.hoomd_snapshot")
        forcefield = import_("foyer.forcefield")
        ff = forcefield.Forcefield(name='oplsaa')
        structure = ff.apply(ethane)
        snap,_ = hoomd_snapshot.to_hoomdsnapshot(structure)

        assert snap.particles.N == 8
        assert snap.bonds.N == 7
        assert snap.angles.N == 12
        assert snap.dihedrals.N == 9
        assert snap.pairs.N == 9

    def test_bad_input_to_hoomdsimulation(self):
        hoomd_simulation = import_("mbuild.formats.hoomd_simulation")
        with pytest.raises(ValueError):
            hoomd_simulation.create_hoomd_simulation('fake_object')

    def test_compound_to_hoomdsimulation(self, ethane):
        hoomd_simulation = import_("mbuild.formats.hoomd_simulation")
        with pytest.raises(ValueError):
            hoomd_simulation.create_hoomd_simulation(ethane)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_structure_to_hoomdsimulation(self, ethane):
        forcefield = import_("foyer.forcefield")
        hoomd = import_("hoomd")
        hoomd_simulation = import_("mbuild.formats.hoomd_simulation")
        ff = forcefield.Forcefield(name='oplsaa')
        structure = ff.apply(ethane)
        hoomd_simulation.create_hoomd_simulation(structure)

        sim_forces = hoomd.context.current.forces
        pair_force = import_("hoomd.md.pair")
        charge_force = import_("hoomd.md.charge")
        special_pair_force = import_("hoomd.md.special_pair")
        bond_force = import_("hoomd.md.bond")
        angle_force = import_("hoomd.md.angle")
        dihedral_force = import_("hoomd.md.dihedral")

        assert isinstance(sim_forces[0], pair_force.lj)
        assert isinstance(sim_forces[1], charge_force.pppm)
        assert isinstance(sim_forces[2], pair_force.ewald)
        assert isinstance(sim_forces[3], special_pair_force.lj)
        assert isinstance(sim_forces[4], special_pair_force.coulomb)
        assert isinstance(sim_forces[5], bond_force.harmonic)
        assert isinstance(sim_forces[6], angle_force.harmonic)
        assert isinstance(sim_forces[7], dihedral_force.opls)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_lj_to_hoomdsimulation(self):
        hoomd = import_("hoomd")
        hoomd_simulation = import_("mbuild.formats.hoomd_simulation")
        forcefield = import_("foyer.forcefield")
        box = mb.Compound()
        box.add(mb.Compound(name='Ar', pos=[1,1,1]))
        box.add(mb.Compound(name='Ar', pos=[1,1,1]))
        ff = forcefield.Forcefield(forcefield_files=get_fn('lj.xml'))
        structure = ff.apply(box)
        structure.box = [10, 10, 10, 90, 90, 90]
        hoomd_simulation.create_hoomd_simulation(structure)
        sim_forces = hoomd.context.current.forces
        pair_force = import_("hoomd.md.pair")

        assert isinstance(sim_forces[0], pair_force.lj)


class TestHoomdXML(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.hoomdxml')

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.hoomdxml', forcefield_name='oplsaa')

    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        ethane.save(filename='ethane-box.hoomdxml', box=box)

    def test_save_triclinic_box_(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]), angles=[60, 70, 80])
        ethane.save(filename='triclinic-box.hoomdxml', box=box)

    def test_rigid(self, benzene):
        n_benzenes = 10
        benzene.name = 'Benzene'
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')
        filled.save(filename='benzene.hoomdxml')

        xml_file = xml.etree.ElementTree.parse('benzene.hoomdxml').getroot()
        body_text = xml_file[0].find('body').text
        rigid_bodies = [int(body) for body in body_text.split('\n') if body]
        for body_id in range(10):
            assert rigid_bodies.count(body_id) == 6
        assert rigid_bodies.count(-1) == n_benzenes * 6

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_number_in_each_section(self, box_of_benzenes):
        box_of_benzenes.save(filename='benzene.hoomdxml', forcefield_name='oplsaa')
        xml_file = xml.etree.ElementTree.parse('benzene.hoomdxml').getroot()
        for attribute in ['position', 'type', 'mass', 'charge']:
            body_text = xml_file[0].find(attribute).text
            list_of_things = [x for x in body_text.split('\n') if x]
            assert len(list_of_things) == 12*10
        for attribute, number in [('bond', 12), ('angle', 18), ('dihedral', 24)]:
            body_text = xml_file[0].find(attribute).text
            list_of_things = [x for x in body_text.split('\n') if x]
            assert len(list_of_things) == number*10

    def test_box_dimensions(self, benzene):
        n_benzenes = 10
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.save(filename='benzene.hoomdxml')
        for atom in mb.load('benzene.hoomdxml'):
            assert atom.pos.max() < 20
            assert atom.pos.min() > -20

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_auto_scale_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.hoomdxml', forcefield_name='oplsaa', auto_scale=True)
        xml_file = xml.etree.ElementTree.parse('ethane-opls.hoomdxml').getroot()
        masses = xml_file[0].find('mass').text.splitlines()
        # We use 1 and 5 since the first element of masses is empty
        assert masses[1] == "1.0"
        assert masses[5] == "1.0"
        pair_coeffs = [_.split("\t") for _ in xml_file[0].find('pair_coeffs').text.splitlines()]
        # The first element is empty, the next element should be ['opls_135', '1.0000', '1.0000']
        assert pair_coeffs[1][1] == "1.0000"
        assert pair_coeffs[1][2] == "1.0000"
