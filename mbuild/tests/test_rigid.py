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

    def test_inherit_rigid(self):
        benzene = mb.Compound()
        rigid_benzene = mb.Compound(rigid=True)

        ch = mb.load(get_fn('ch.mol2'))
        mb.translate(ch, -ch[0].pos)        
        ch.add(mb.Port(anchor=ch[0]), 'a')
        mb.translate(ch['a'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['a'], 120.0 * (np.pi/180.0))

        ch.add(mb.Port(anchor=ch[0]), 'b')
        mb.translate(ch['b'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['b'], -120.0 * (np.pi/180.0))

        ch_clone = mb.clone(ch)
        benzene.add(ch)
        rigid_benzene.add(ch_clone)
        current1 = ch
        current2 = ch_clone

        for _ in range(5):
            ch_new = mb.clone(ch)
            ch_new_clone = mb.clone(ch_clone)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current1['b'])
            mb.force_overlap(move_this=ch_new_clone,
                             from_positions=ch_new_clone['a'],
                             to_positions=current2['b'])
            current1 = ch_new
            current2 = ch_new_clone
            benzene.add(ch_new)
            rigid_benzene.add(ch_new_clone)

        carbons1 = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons1[0],carbons1[-1]))
        carbons2 = [p for p in rigid_benzene.particles_by_name('C')]
        rigid_benzene.add_bond((carbons2[0],carbons2[-1]))

        benzene.inherit_rigid()
        rigid_benzene.inherit_rigid()

        assert benzene.rigid is False
        assert rigid_benzene.rigid is 0
        assert not [p for p in benzene.rigid_particles()]
        assert [p for p in rigid_benzene.rigid_particles()]
        assert not any(v is 0 for v in [p for p in benzene.rigid_ids()])
        assert all(v is 0 for v in [p for p in rigid_benzene.rigid_ids()])
