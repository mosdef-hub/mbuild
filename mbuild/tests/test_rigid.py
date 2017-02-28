import pytest
import numpy as np

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn


class TestRigid(BaseTest):

    def test_load_rigid(self, rigid_benzene):
        assert rigid_benzene.rigid is True
        assert rigid_benzene.rigid_id is 0
        assert rigid_benzene.max_rigid_id is 0
        assert len(list(rigid_benzene.rigid_particles())) == 12
        assert [p for p in rigid_benzene.rigid_ids()].count(0) == 12

    def test_load_nonrigid(self, benzene):
        assert benzene.rigid is False
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is None
        assert len(list(benzene.rigid_particles())) == 0
        assert [p for p in benzene.rigid_ids()].count(None) == 12

    def test_rigid_from_parts(self, rigid_ch):
        benzene = mb.Compound()
        benzene.add(rigid_ch)
        current = rigid_ch

        for _ in range(5):
            ch_new = mb.clone(rigid_ch)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current['b'])
            current = ch_new
            benzene.add(ch_new, increment_rigid=False)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        assert benzene.rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert len(list(benzene.rigid_particles())) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12

    def test_rigid_from_parts2(self, rigid_ch):
        benzene = mb.Compound(rigid=True)
        benzene.add(rigid_ch, increment_rigid=False)
        current = rigid_ch

        for _ in range(5):
            ch_new = mb.clone(rigid_ch)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current['b'])
            current = ch_new
            benzene.add(ch_new, increment_rigid=False)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        assert benzene.rigid is True
        assert benzene.rigid_id is 0
        assert benzene.max_rigid_id is 0
        assert len(list(benzene.rigid_particles())) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12

    def test_nonrigid_from_parts(self, benzene_from_parts):
        assert benzene_from_parts.rigid is False
        assert benzene_from_parts.rigid_id is None
        assert benzene_from_parts.max_rigid_id is None
        assert len(list(benzene_from_parts.rigid_particles())) == 0
        assert [p for p in benzene_from_parts.rigid_ids()].count(None) == 12

    def test_set_rigid(self, benzene):
        benzene.set_rigid()

        assert benzene.rigid_id is 0
        assert benzene.max_rigid_id is 0
        assert len(list(benzene.rigid_particles())) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12

    def test_set_rigid_by_name(self, benzene):
        benzene.set_rigid(name='C')

        assert benzene.rigid_id is 0
        assert benzene.max_rigid_id is 0
        assert len(list(benzene.rigid_particles())) == 6
        assert [p for p in benzene.rigid_ids()].count(0) == 6
        assert [p for p in benzene.rigid_ids()].count(False) == 6

    def test_increment_rigid_id(self, rigid_benzene):
        compound = mb.Compound()
        rigid_benzene2 = mb.clone(rigid_benzene)
        compound.add(rigid_benzene)
        compound.add(rigid_benzene2)

        assert rigid_benzene.rigid_id is 0
        assert rigid_benzene2.rigid_id is 1
        assert compound.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 24
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12

    def test_increment_rigid_id_partial(self, benzene):
        compound = mb.Compound()
        benzene.set_rigid(name='C')
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2, increment_rigid=True)

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 12
        assert [p for p in compound.rigid_ids()].count(0) == 6
        assert [p for p in compound.rigid_ids()].count(1) == 6
        assert [p for p in compound.rigid_ids()].count(None) == 12

    def test_turn_into_rigid_mix(self, benzene):
        compound = mb.Compound()
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        benzene.set_rigid()

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is None
        assert compound.max_rigid_id is 0
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(None) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(None) == 12

    def test_turn_into_rigid_individual(self, benzene):
        compound = mb.Compound()
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        benzene.set_rigid()
        benzene2.set_rigid()

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid_id is 1
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(1) == 12

    def test_turn_into_rigid_multiple(self, benzene):
        compound = mb.Compound()
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.set_rigid()

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 0
        assert compound.max_rigid_id is 0
        assert [p for p in compound.rigid_ids()].count(0) == 24
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(0) == 12

    def test_create_rigid_bodies(self, benzene):
        compound = mb.Compound()
        benzene.name = 'Benzene'
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.create_rigid_bodies(name='Benzene')

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid_id is 1
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(1) == 12

    def test_create_rigid_bodies_child_name(self, benzene):
        compound = mb.Compound()
        benzene.name = 'Benzene'
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.create_rigid_bodies(name='Benzene', particle_name='C')

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid_id is 1
        assert [p for p in compound.rigid_ids()].count(0) == 6
        assert [p for p in compound.rigid_ids()].count(1) == 6
        assert [p for p in compound.rigid_ids()].count(None) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 6
        assert [p for p in benzene2.rigid_ids()].count(1) == 6

    def test_fill_box_rigid(self, rigid_benzene):
        n_benzenes = 10
        filled = mb.fill_box(rigid_benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])

        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * rigid_benzene.n_particles

    def test_fill_box_semi_rigid(self, benzene):
        n_benzenes = 10
        benzene.set_rigid(name='C')
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])

        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 6

    def test_set_create_rigid_bodies_after_fill(self, benzene):
        n_benzenes = 10
        benzene.name = 'Benzene'
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene')

        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * benzene.n_particles

    def test_set_create_rigid_bodies_duplicate_warn(self, rigid_benzene):
        with pytest.warns(UserWarning):
            n_benzenes = 10
            rigid_benzene.name= 'Benzene'
            filled = mb.fill_box(rigid_benzene,
                                 n_compounds=n_benzenes,
                                 box=[0, 0, 0, 4, 4, 4])
            filled.create_rigid_bodies(name='Benzene')

    def test_create_semi_rigid_bodies_after_fill(self, benzene):
        n_benzenes = 10
        benzene.name = 'Benzene'
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene', particle_name='C')

        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 6

    def test_create_semi_rigid_bodies_hierarchy(self, benzene_from_parts):
        n_benzenes = 10
        benzene_from_parts.name = 'Benzene'
        filled = mb.fill_box(benzene_from_parts,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene', particle_name='C')

        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 6

    def test_create_semi_rigid_bodies_hierarchy2(self, benzene_from_parts):
        n_benzenes = 10
        benzene_from_parts.name = 'Benzene'
        filled = mb.fill_box(benzene_from_parts,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled2 = mb.clone(filled)
        compound = mb.Compound(subcompounds=[filled, filled2])
        compound.create_rigid_bodies(name='Benzene', particle_name='C')

        assert compound.max_rigid_id == (n_benzenes*2) - 1
        assert len(list(compound.rigid_particles())) == n_benzenes * 2 * 6

    def test_create_semi_rigid_bodies_hierarchy3(self, benzene_from_parts):
        n_benzenes = 10
        benzene_from_parts.name = 'Benzene'
        filled = mb.fill_box(benzene_from_parts,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene', particle_name='C')
        filled2 = mb.clone(filled)
        filled.add(filled2)

        assert filled.max_rigid_id == (n_benzenes*2) - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 2 * 6
        assert [p for p in filled.rigid_ids()].count(1) == 6

    def test_create_semi_rigid_bodies_hierarchy4(self, benzene_from_parts):
        n_benzenes = 10
        benzene_from_parts.name = 'Benzene'
        filled = mb.fill_box(benzene_from_parts,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene', particle_name='C')
        filled2 = mb.clone(filled)
        filled.add(filled2, increment_rigid=False)

        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 2 * 6
        assert [p for p in filled.rigid_ids()].count(1) == 12

    def test_delete_body(self, rigid_benzene):
        n_benzenes = 10
        filled = mb.fill_box(rigid_benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.remove(filled.children[0])

        assert filled.max_rigid_id == n_benzenes - 2
        assert len(list(filled.rigid_particles())) == (n_benzenes - 1) * rigid_benzene.n_particles

    def test_delete_body_semi_rigid(self, benzene):
        n_benzenes = 10
        benzene.name = 'Benzene'
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene', particle_name='C')
        filled.remove(filled.children[0])

        assert filled.max_rigid_id == n_benzenes - 2
        assert len(list(filled.rigid_particles())) == (n_benzenes - 1) * 6

    def test_rigid_with_subcompounds1(self, benzene):
        compound = mb.Compound(subcompounds=benzene, rigid=True)
        assert compound.rigid is True
        assert benzene.rigid_id is 0
        assert compound.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 12
        assert all(v is 0 for v in [p for p in compound.rigid_ids()])

    def test_rigid_with_subcompounds2(self, benzene):
        benzene2 = mb.clone(benzene)
        compound = mb.Compound(subcompounds=[benzene, benzene2], rigid=True)

        assert compound.rigid is True
        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 0
        assert len(list(compound.rigid_particles())) == 24
        assert [p for p in compound.rigid_ids()].count(0) == 24
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(0) == 12

    def test_rigid_with_subcompounds3(self, rigid_benzene):
        compound = mb.Compound(subcompounds=rigid_benzene)
        assert compound.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 12
        assert all(v is 0 for v in [p for p in compound.rigid_ids()])

    def test_rigid_with_subcompounds4(self, benzene):
        benzene.set_rigid(name='C')
        compound = mb.Compound(subcompounds=benzene)
        assert compound.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 6

    def test_rigid_with_subcompounds5(self, benzene):
        benzene.set_rigid(name='C')
        benzene2 = mb.clone(benzene)
        compound = mb.Compound(subcompounds=[benzene, benzene2])

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 12
        assert [p for p in compound.rigid_ids()].count(0) == 6
        assert [p for p in compound.rigid_ids()].count(1) == 6
