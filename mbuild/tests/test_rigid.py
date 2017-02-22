import pytest
import numpy as np
import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn


class TestRigid(BaseTest):

    def test_load_rigid(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        assert benzene.rigid is True
        assert benzene.rigid_id is 0
        assert benzene.max_rigid() is 0
        assert len([p for p in benzene.rigid_particles()]) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12

    def test_load_nonrigid(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        assert benzene.rigid is False
        assert benzene.rigid_id is None
        assert benzene.max_rigid() is None
        assert len([p for p in benzene.rigid_particles()]) == 0
        assert [p for p in benzene.rigid_ids()].count(None) == 12

    def test_rigid_from_parts(self):
        benzene = mb.Compound()

        ch = mb.load(get_fn('ch.mol2'), rigid=True)
        mb.translate(ch, -ch[0].pos)        
        ch.add(mb.Port(anchor=ch[0]), 'a')
        mb.translate(ch['a'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['a'], 120.0 * (np.pi/180.0))

        ch.add(mb.Port(anchor=ch[0]), 'b')
        mb.translate(ch['b'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['b'], -120.0 * (np.pi/180.0))

        benzene.add(ch)
        current = ch

        for _ in range(5):
            ch_new = mb.clone(ch)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current['b'])
            current = ch_new
            benzene.add(ch_new, increment_rigid=False)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        assert benzene.rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid() is 0
        assert len([p for p in benzene.rigid_particles()]) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12

    def test_rigid_from_parts2(self):
        benzene = mb.Compound(rigid=True)

        ch = mb.load(get_fn('ch.mol2'), rigid=True)
        mb.translate(ch, -ch[0].pos)        
        ch.add(mb.Port(anchor=ch[0]), 'a')
        mb.translate(ch['a'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['a'], 120.0 * (np.pi/180.0))

        ch.add(mb.Port(anchor=ch[0]), 'b')
        mb.translate(ch['b'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['b'], -120.0 * (np.pi/180.0))

        benzene.add(ch)
        current = ch

        for _ in range(5):
            ch_new = mb.clone(ch)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current['b'])
            current = ch_new
            benzene.add(ch_new, increment_rigid=False)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        assert benzene.rigid is True
        assert benzene.rigid_id is 0
        assert benzene.max_rigid() is 0
        assert len([p for p in benzene.rigid_particles()]) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12

    def test_nonrigid_from_parts(self):
        benzene = mb.Compound()

        ch = mb.load(get_fn('ch.mol2'))
        mb.translate(ch, -ch[0].pos)        
        ch.add(mb.Port(anchor=ch[0]), 'a')
        mb.translate(ch['a'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['a'], 120.0 * (np.pi/180.0))

        ch.add(mb.Port(anchor=ch[0]), 'b')
        mb.translate(ch['b'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['b'], -120.0 * (np.pi/180.0))

        benzene.add(ch)
        current = ch

        for _ in range(5):
            ch_new = mb.clone(ch)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current['b'])
            current = ch_new
            benzene.add(ch_new)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        assert benzene.rigid is False
        assert benzene.rigid_id is None
        assert benzene.max_rigid() is None
        assert len([p for p in benzene.rigid_particles()]) == 0
        assert [p for p in benzene.rigid_ids()].count(None) == 12

    def test_set_rigid(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.set_rigid()

        assert benzene.rigid_id is 0
        assert benzene.max_rigid() is 0
        assert len([p for p in benzene.rigid_particles()]) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12

    def test_set_rigid_by_name(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.set_rigid(name='C')

        assert benzene.rigid_id is 0
        assert benzene.max_rigid() is 0
        assert len([p for p in benzene.rigid_particles()]) == 6
        assert [p for p in benzene.rigid_ids()].count(0) == 6
        assert [p for p in benzene.rigid_ids()].count(False) == 6

    def test_increment_rigid_id(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid() is 1
        assert len([p for p in compound.rigid_particles()]) == 24
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12

    def test_increment_rigid_id_partial(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.set_rigid(name='C')
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2, increment_rigid=True)

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid() is 1
        assert len([p for p in compound.rigid_particles()]) == 12
        assert [p for p in compound.rigid_ids()].count(0) == 6
        assert [p for p in compound.rigid_ids()].count(1) == 6
        assert [p for p in compound.rigid_ids()].count(None) == 12

    def test_turn_into_rigid_mix(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        benzene.set_rigid()

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is None
        assert compound.max_rigid() is 0
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(None) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(None) == 12

    def test_turn_into_rigid_individual(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        benzene.set_rigid()
        benzene2.set_rigid()

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid() is 1
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(1) == 12

    def test_turn_into_rigid_multiple(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.set_rigid()

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 0
        assert compound.max_rigid() is 0
        assert [p for p in compound.rigid_ids()].count(0) == 24
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(0) == 12

    def test_create_rigid_bodies(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.name = 'Benzene'
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.create_rigid_bodies(name='Benzene')

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid() is 1
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(1) == 12

    def test_create_rigid_bodies_child_name(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.name = 'Benzene'
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.create_rigid_bodies(name='Benzene', particle_name='C')

        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert compound.max_rigid() is 1
        assert [p for p in compound.rigid_ids()].count(0) == 6
        assert [p for p in compound.rigid_ids()].count(1) == 6
        assert [p for p in compound.rigid_ids()].count(None) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 6
        assert [p for p in benzene2.rigid_ids()].count(1) == 6

    def test_fill_box_rigid(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        n_benzenes = 10
        filled = mb.fill_box(benzene, n_compounds=n_benzenes, box=[0, 0, 0, 4, 4, 4])

        assert filled.max_rigid() == n_benzenes - 1
        assert len([p for p in filled.rigid_particles()]) == n_benzenes * benzene.n_particles

    def test_fill_box_semi_rigid(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        n_benzenes = 10
        benzene.set_rigid(name='C')
        filled = mb.fill_box(benzene, n_compounds=n_benzenes, box=[0, 0, 0, 4, 4, 4])

        assert filled.max_rigid() == n_benzenes - 1
        assert len([p for p in filled.rigid_particles()]) == n_benzenes * 6

    def test_set_create_rigid_bodies_after_fill(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        n_benzenes = 10
        benzene.name= 'Benzene'
        filled = mb.fill_box(benzene, n_compounds=n_benzenes, box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene')

        assert filled.max_rigid() == n_benzenes - 1
        assert len([p for p in filled.rigid_particles()]) == n_benzenes * benzene.n_particles

    def test_create_semi_rigid_bodies_after_fill(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        n_benzenes = 10
        benzene.name = 'Benzene'
        filled = mb.fill_box(benzene, n_compounds=10, box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene', particle_name='C')

        assert filled.max_rigid() == n_benzenes - 1
        assert len([p for p in filled.rigid_particles()]) == n_benzenes * 6

    def test_delete_body(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        n_benzenes = 10
        filled = mb.fill_box(benzene, n_compounds=n_benzenes, box=[0, 0, 0, 4, 4, 4])
        filled.remove(filled.children[0])

        assert filled.max_rigid() == n_benzenes - 2
        assert len([p for p in filled.rigid_particles()]) == (n_benzenes - 1) * benzene.n_particles

    def test_delete_body_semi_rigid(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        n_benzenes = 10
        benzene.name = 'Benzene'
        filled = mb.fill_box(benzene, n_compounds=10, box=[0, 0, 0, 4, 4, 4])
        filled.create_rigid_bodies(name='Benzene', particle_name='C')
        filled.remove(filled.children[0])

        assert filled.max_rigid() == n_benzenes - 2
        assert len([p for p in filled.rigid_particles()]) == (n_benzenes - 1) * 6

    def test_rigid_with_subcompounds1(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        compound = mb.Compound(subcompounds=benzene, rigid=True)
        assert compound.rigid is True
        assert benzene.rigid_id is 0
        assert compound.max_rigid() is 0
        assert len([p for p in compound.rigid_particles()]) == 12
        assert all(v is 0 for v in [p for p in compound.rigid_ids()])

    def test_rigid_with_subcompounds2(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        compound = mb.Compound(subcompounds=[benzene, benzene2])
        assert compound.rigid is True
        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert len([p for p in compound.rigid_particles()]) == 24
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(1) == 12

    def test_rigid_with_subcompounds3(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        benzene3 = mb.clone(benzene)
        compound = mb.Compound([benzene, [benzene2, benzene3]])
        assert compound.rigid is True
        assert benzene.rigid_id is 0
        assert benzene2.rigid_id is 1
        assert benzene3.rigid_id is 2
        assert len([p for p in compound.rigid_particles()]) == 36
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12
        assert [p for p in compound.rigid_ids()].count(2) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(1) == 12
        assert [p for p in benzene3.rigid_ids()].count(2) == 12

    def test_rigid_with_subcompounds4(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.set_rigid(name='C')
        compound = mb.Compound(benzene, rigid=True)

    def test_rigid_with_subcompounds5(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.set_rigid(name='C')
        benzene2 = mb.clone(benzene)
        compound = mb.Compound([benzene, benzene2])

    def test_rigid_with_subcompounds6(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        compound = mb.Compound(benzene, rigid=True)

    def test_rigid_with_subcompounds7(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        benzene2 = mb.clone(benzene)
        compound = mb.Compound([benzene, benzene2], rigid=True)

    def test_rigid_with_subcompounds8(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        compound = mb.Compound(benzene)

    def test_rigid_with_subcompounds9(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        benzene2 = mb.clone(benzene)
        compound = mb.Compound([benzene, benzene2])
