import pytest
import numpy as np
import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn


class TestRigid(BaseTest):

    def test_load_rigid(self):
        rigid_benzene = mb.load(get_fn('benzene.mol2'), rigid=True)

        assert rigid_benzene.rigid is 0

    def test_load_nonrigid(self):
        benzene = mb.load(get_fn('benzene.mol2'))

        assert benzene.rigid is False

    def test_rigid_from_parts(self):
        rigid_benzene = mb.Compound(rigid=True)

        ch = mb.load(get_fn('ch.mol2'))
        mb.translate(ch, -ch[0].pos)        
        ch.add(mb.Port(anchor=ch[0]), 'a')
        mb.translate(ch['a'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['a'], 120.0 * (np.pi/180.0))

        ch.add(mb.Port(anchor=ch[0]), 'b')
        mb.translate(ch['b'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['b'], -120.0 * (np.pi/180.0))

        rigid_benzene.add(ch)
        current = ch

        for _ in range(5):
            ch_new = mb.clone(ch)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current1['b'])
            current = ch_new
            rigid_benzene.add(ch_new)

        carbons = [p for p in rigid_benzene.particles_by_name('C')]
        rigid_benzene.add_bond((carbons2[0],carbons2[-1]))

        rigid_benzene.set_rigid()

        assert rigid_benzene.rigid is 0
        assert [p for p in rigid_benzene.rigid_particles()]
        assert all(v is 0 for v in [p for p in rigid_benzene.rigid_ids()])

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
                             to_positions=current1['b'])
            current = ch_new
            benzene.add(ch_new)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons2[0],carbons2[-1]))

        benzene.set_rigid()

        assert benzene.rigid is False
        assert not [p for p in benzene.rigid_particles()]
        assert all(v is False for v in [p for p in benzene.rigid_ids()])

    def test_set_rigid_by_name(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.set_rigid(name='C')

        assert len([p for p in benzene.rigid_particles()]) == 6
        assert benzene.rigid is 0
        assert [p for p in benzene.rigid_ids()].count(0) == 6
        assert [p for p in benzene.rigid_ids()].count(False) == 6

    def test_increment_rigid_id(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)

        assert benzene.rigid is 0
        assert benzene2.rigid is 1
        assert compound.rigid is True

        rigid1 = np.zeros(len(benzene.n_particles))
        rigid2 = np.ones(len(benzene2.n_particles))
        assert np.array_equal([p for p in compound.rigid_ids()],np.hstack((rigid1,rigid2)))

    def test_increment_rigid_id_partial(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.set_rigid(name='C')
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)

        assert benzene.rigid is 0
        assert benzene2.rigid is 1
        assert compound.rigid is True
        assert [p for p in compound.rigid_ids()].count(0) == 6
        assert [p for p in compound.rigid_ids()].count(1) == 6
        assert [p for p in compound.rigid_ids()].count(False) == 12

    def test_turn_into_rigid(self):
        rigid_benzene = mb.load(get_fn('benzene.mol2'))
        rigid_benzene.rigid = True
        assert all(v is 0 for v in [p for p in rigid_benzene.rigid_ids()])

    def test_turn_into_rigid_mix(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        benzene.rigid = True

        assert benzene.rigid is 0
        assert benzene2.rigid is False
        assert compound.rigid is True
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(False) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(False) == 12

    def test_turn_into_rigid_individual(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        benzene.rigid = True
        benzene2.rigid = True

        assert benzene.rigid is 0
        assert benzene2.rigid is 1
        assert compound.rigid is True
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
        compound.rigid = True

        assert benzene.rigid is 0
        assert benzene2.rigid is 0
        assert compound.rigid is True
        assert [p for p in compound.rigid_ids()].count(0) == 24
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(0) == 12

    def test_set_rigid_parent_name(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.set_rigid(parent_name='Benzene')

        assert benzene.rigid is 0
        assert benzene2.rigid is 1
        assert compound.rigid is True
        assert [p for p in compound.rigid_ids()].count(0) == 12
        assert [p for p in compound.rigid_ids()].count(1) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 12
        assert [p for p in benzene2.rigid_ids()].count(1) == 12

    def test_set_rigid_particle_name_parent_name(self):
        compound = mb.Compound()
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.set_rigid(particle_name='C', parent_name='Benzene')

        assert benzene.rigid is 0
        assert benzene2.rigid is 1
        assert compound.rigid is True
        assert [p for p in compound.rigid_ids()].count(0) == 6
        assert [p for p in compound.rigid_ids()].count(1) == 6
        assert [p for p in compound.rigid_ids()].count(False) == 12
        assert [p for p in benzene.rigid_ids()].count(0) == 6
        assert [p for p in benzene.rigid_ids()].count(False) == 6
        assert [p for p in benzene2.rigid_ids()].count(1) == 6
        assert [p for p in benzene2.rigid_ids()].count(False) == 6

    def test_fill_box_rigid(self):
        benzene = mb.load(get_fn('benzene.mol2'), rigid=True)
        filled = mb.fill_box(benzene, n_compounds=10, box=[0, 0, 0, 4, 4, 4])

    def test_fill_box_semi_rigid(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        benzene.set_rigid(particle_name='C')
        filled = mb.fill_box(benzene, n_compounds=10, box=[0, 0, 0, 4, 4, 4])

    def test_set_rigid_after_fill(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        filled = mb.fill_box(benzene, n_compounds=10, box=[0, 0, 0, 4, 4, 4])
        filled.set_rigid(parent_name='Benzene')

    def test_set_rigid_after_fill_partial(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        filled = mb.fill_box(benzene, n_compounds=10, box=[0, 0, 0, 4, 4, 4])
        filled.set_rigid(particle_name='C', parent_name='Benzene')

    def test_delete_body(self):
        benzene = mb.load(get_fn('benzene.mol2'))
        filled = mb.fill_box(benzene, n_compounds=10, box=[0, 0, 0, 4, 4, 4])
        filled.set_rigid(parent_name='Benzene')
        filled.remove(filled.children[0])
