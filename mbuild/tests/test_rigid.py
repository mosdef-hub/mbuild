import pytest
import numpy as np

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn


class TestRigid(BaseTest):

    def test_load_rigid(self, rigid_benzene):
        assert rigid_benzene.contains_rigid is True
        assert rigid_benzene[0].contains_rigid is False
        assert rigid_benzene.rigid_id is None
        assert rigid_benzene.max_rigid_id is 0
        assert len(list(rigid_benzene.rigid_particles(rigid_id=0))) == 12

    def test_load_nonrigid(self, benzene):
        assert benzene.contains_rigid is False
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is None
        assert len(list(benzene.rigid_particles())) == 0

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
            benzene.add(ch_new, reset_rigid_ids=False)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene.children[0].contains_rigid == True
        assert benzene.children[0].rigid_id is None
        assert len(list(benzene.rigid_particles(rigid_id=0))) == 12

    def test_rigid_from_parts2(self, rigid_ch):
        benzene = mb.Compound()
        benzene.add(rigid_ch, reset_rigid_ids=False)
        current = rigid_ch

        for _ in range(5):
            ch_new = mb.clone(rigid_ch)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current['b'])
            current = ch_new
            benzene.add(ch_new, reset_rigid_ids=False)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene.children[0].contains_rigid == True
        assert benzene.children[0].rigid_id is None
        assert len(list(benzene.rigid_particles(rigid_id=0))) == 12

    def test_rigid_from_parts3(self, rigid_ch):
        benzene = mb.Compound()
        benzene.add(rigid_ch)
        current = rigid_ch

        for _ in range(5):
            ch_new = mb.clone(rigid_ch)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current['b'])
            current = ch_new
            benzene.add(ch_new)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 5
        assert benzene.children[0].contains_rigid == True
        assert benzene.children[0].rigid_id is None
        for rigid_id in range(6):
            assert len(list(benzene.rigid_particles(rigid_id=rigid_id))) == 2

    def test_nonrigid_from_parts(self, benzene_from_parts):
        assert benzene_from_parts.contains_rigid is False
        assert benzene_from_parts.rigid_id is None
        assert benzene_from_parts.max_rigid_id is None
        assert len(list(benzene_from_parts.rigid_particles())) == 0

    def test_label_rigid_bodies_single_partial(self, benzene):
        benzene.label_rigid_bodies(rigid_particles='C')

        assert benzene.contains_rigid == True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert len(list(benzene.rigid_particles())) == 6
        assert len(list(benzene.rigid_particles(rigid_id=0))) == 6

    def test_save_non_sequential_rigid_ids(self, benzene):
        n_benzenes = 10
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.label_rigid_bodies(discrete_bodies='Benzene')
        filled.children[0]._increment_rigid_ids(increment=3)
        with pytest.warns(UserWarning):
            filled.save('benzene-box.hoomdxml')

    def test_increment_rigid_id(self, rigid_benzene):
        compound = mb.Compound()
        rigid_benzene2 = mb.clone(rigid_benzene)
        compound.add(rigid_benzene)
        compound.add(rigid_benzene2)

        assert rigid_benzene.contains_rigid is True
        assert rigid_benzene.rigid_id is None
        assert rigid_benzene.max_rigid_id is 0
        assert rigid_benzene2.contains_rigid is True
        assert rigid_benzene2.rigid_id is None
        assert rigid_benzene2.max_rigid_id is 1
        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 1
        assert len(list(compound.rigid_particles(rigid_id=0))) == 12
        assert len(list(compound.rigid_particles(rigid_id=1))) == 12

    def test_increment_rigid_id_partial(self, benzene):
        compound = mb.Compound()
        benzene.label_rigid_bodies(rigid_particles='C')
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene2.contains_rigid is True
        assert benzene2.rigid_id is None
        assert benzene2.max_rigid_id is 1
        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 12
        assert len(list(compound.rigid_particles(rigid_id=0))) == 6
        assert len(list(compound.rigid_particles(rigid_id=1))) == 6

    def test_turn_into_rigid_after_add(self, benzene):
        compound = mb.Compound()
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        benzene.label_rigid_bodies()

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene2.contains_rigid is False
        assert benzene2.rigid_id is None
        assert benzene2.max_rigid_id is None
        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 12
        assert len(list(compound.rigid_particles(rigid_id=0))) == 12
        assert len(list(benzene.rigid_particles(rigid_id=0))) == 12
        assert len(list(benzene2.rigid_particles())) == 0

    def test_turn_into_rigid_after_add_multi(self, benzene):
        compound = mb.Compound()
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        benzene.label_rigid_bodies()
        benzene2.label_rigid_bodies()

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene2.contains_rigid is True
        assert benzene2.rigid_id is None
        assert benzene2.max_rigid_id is 1
        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 24
        assert len(list(compound.rigid_particles(rigid_id=0))) == 12
        assert len(list(compound.rigid_particles(rigid_id=1))) == 12
        assert len(list(benzene.rigid_particles(rigid_id=0))) == 12
        assert len(list(benzene2.rigid_particles(rigid_id=1))) == 12

    def test_turn_into_rigid_after_add_parent(self, benzene):
        compound = mb.Compound()
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.label_rigid_bodies()

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene2.contains_rigid is True
        assert benzene2.rigid_id is None
        assert benzene2.max_rigid_id is 0
        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 24
        assert len(list(compound.rigid_particles(rigid_id=0))) == 24
        assert len(list(benzene.rigid_particles(rigid_id=0))) == 12
        assert len(list(benzene2.rigid_particles(rigid_id=0))) == 12

    def test_label_rigid_bodies_multi(self, benzene):
        compound = mb.Compound()
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.label_rigid_bodies(discrete_bodies='Benzene')

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene2.contains_rigid is True
        assert benzene2.rigid_id is None
        assert benzene2.max_rigid_id is 1
        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 24
        assert len(list(compound.rigid_particles(rigid_id=0))) == 12
        assert len(list(compound.rigid_particles(rigid_id=1))) == 12
        assert len(list(benzene.rigid_particles(rigid_id=0))) == 12
        assert len(list(benzene2.rigid_particles(rigid_id=1))) == 12

    def test_label_rigid_bodies_multi_partial(self, benzene):
        compound = mb.Compound()
        benzene2 = mb.clone(benzene)
        compound.add(benzene)
        compound.add(benzene2)
        compound.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')

        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene2.contains_rigid is True
        assert benzene2.rigid_id is None
        assert benzene2.max_rigid_id is 1
        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 12
        assert len(list(compound.rigid_particles(rigid_id=0))) == 6
        assert len(list(compound.rigid_particles(rigid_id=1))) == 6
        assert len(list(benzene.rigid_particles(rigid_id=0))) == 6
        assert len(list(benzene2.rigid_particles(rigid_id=1))) == 6

    def test_fill_box_rigid(self, rigid_benzene):
        n_benzenes = 10
        filled = mb.fill_box(rigid_benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])

        assert filled.contains_rigid is True
        assert filled.rigid_id is None
        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * rigid_benzene.n_particles

    def test_fill_box_semi_rigid(self, benzene):
        n_benzenes = 10
        benzene.label_rigid_bodies(rigid_particles='C')
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])

        assert filled.contains_rigid is True
        assert filled.rigid_id is None
        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 6

    def test_label_rigid_bodies_after_fill(self, benzene):
        n_benzenes = 10
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.label_rigid_bodies(discrete_bodies='Benzene')

        assert filled.contains_rigid is True
        assert filled.rigid_id is None
        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * benzene.n_particles

    def test_label_rigid_bodies_list(self, benzene):
        n_benzenes = 10
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.children[0].name = 'Benzene0'
        filled.children[1].name = 'Benzene1'
        filled.label_rigid_bodies(discrete_bodies=['Benzene0', 'Benzene1'])

        assert filled.contains_rigid is True
        assert filled.rigid_id is None
        assert filled.max_rigid_id == 1
        assert len(list(filled.rigid_particles())) == 2 * benzene.n_particles

    def test_label_rigid_bodies_list_particle_list(self, benzene):
        n_benzenes = 10
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.children[0].name = 'Benzene0'
        filled.children[1].name = 'Benzene1'
        filled.label_rigid_bodies(discrete_bodies=['Benzene0', 'Benzene1'],
                                  rigid_particles=['C', 'H'])

        assert filled.contains_rigid is True
        assert filled.rigid_id is None
        assert filled.max_rigid_id == 1
        assert len(list(filled.rigid_particles())) == 2 * benzene.n_particles

    def test_label_rigid_bodies_duplicate_warn(self, rigid_benzene):
        with pytest.warns(UserWarning):
            n_benzenes = 10
            filled = mb.fill_box(rigid_benzene,
                                 n_compounds=n_benzenes,
                                 box=[0, 0, 0, 4, 4, 4])
            filled.label_rigid_bodies(discrete_bodies='Benzene')

    def test_label_semi_rigid_bodies_after_fill(self, benzene):
        n_benzenes = 10
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')

        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 6

    def test_create_semi_rigid_bodies_hierarchy(self, benzene_from_parts):
        n_benzenes = 10
        filled = mb.fill_box(benzene_from_parts,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.name = 'Benzene box'
        filled2 = mb.clone(filled)
        compound = mb.Compound(subcompounds=[filled, filled2])

        compound.label_rigid_bodies(discrete_bodies='Benzene box')
        assert compound.max_rigid_id == 1
        assert filled.max_rigid_id == 0
        assert filled2.max_rigid_id == 1
        assert len(list(compound.rigid_particles())) == n_benzenes * 2 * 12

        compound.unlabel_rigid_bodies()
        compound.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')
        assert compound.max_rigid_id == (n_benzenes*2) - 1
        assert filled.max_rigid_id == n_benzenes - 1
        assert filled2.max_rigid_id == (n_benzenes*2) - 1
        assert len(list(compound.rigid_particles())) == n_benzenes * 2 * 6
        assert len(list(filled.rigid_particles())) == n_benzenes * 6
        assert len(list(filled2.rigid_particles())) == n_benzenes * 6

        compound.unlabel_rigid_bodies()
        compound.label_rigid_bodies(discrete_bodies='CH')
        assert compound.max_rigid_id == (n_benzenes*2*6) - 1
        assert filled.max_rigid_id == (n_benzenes*6) - 1
        assert filled2.max_rigid_id == (n_benzenes*2*6) - 1
        assert len(list(compound.rigid_particles())) == n_benzenes * 2 * 12
        assert len(list(filled.rigid_particles())) == n_benzenes * 12
        assert len(list(filled2.rigid_particles())) == n_benzenes * 12

    def test_create_semi_rigid_bodies_filled_clone(self, benzene_from_parts):
        n_benzenes = 10
        filled = mb.fill_box(benzene_from_parts,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')
        filled2 = mb.clone(filled)
        filled.add(filled2)

        assert filled.max_rigid_id == (n_benzenes*2) - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 2 * 6
        for rigid_id in range(n_benzenes * 2):
            assert len(list(filled.rigid_particles(rigid_id=rigid_id))) == 6

    def test_create_semi_rigid_bodies_filled_no_increment(self, benzene_from_parts):
        n_benzenes = 10
        filled = mb.fill_box(benzene_from_parts,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')
        filled2 = mb.clone(filled)
        filled.add(filled2, reset_rigid_ids=False)

        assert filled.max_rigid_id == n_benzenes - 1
        assert len(list(filled.rigid_particles())) == n_benzenes * 2 * 6
        for rigid_id in range(n_benzenes):
            assert len(list(filled.rigid_particles(rigid_id=rigid_id))) == 12

    def test_delete_body(self, rigid_benzene):
        n_benzenes = 10
        filled = mb.fill_box(rigid_benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.remove(filled.children[0])

        assert filled.max_rigid_id == n_benzenes - 2
        assert len(list(filled.rigid_particles())) == (n_benzenes - 1) * rigid_benzene.n_particles

    def test_delete_body_particle_by_particle(self, rigid_benzene):
        n_benzenes = 10
        filled = mb.fill_box(rigid_benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        for particle in filled.children[0].particles():
            filled.remove(particle)

        assert filled.max_rigid_id == n_benzenes - 2
        assert len(list(filled.rigid_particles())) == (n_benzenes - 1) * rigid_benzene.n_particles

    def test_delete_body_multiple(self, rigid_benzene):
        n_benzenes = 10
        filled = mb.fill_box(rigid_benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.remove([filled.children[0], filled.children[1]])

        assert filled.max_rigid_id == n_benzenes - 3
        assert len(list(filled.rigid_particles())) == (n_benzenes - 2) * rigid_benzene.n_particles

    def test_delete_body_all(self, rigid_benzene):
        n_benzenes = 10
        filled = mb.fill_box(rigid_benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        for i, child in enumerate(filled.children[:-1]):
            filled.remove(child)
            assert filled.max_rigid_id == n_benzenes - 1 - (i + 1)
            assert len(list(filled.rigid_particles())) == (n_benzenes - (i + 1)) * rigid_benzene.n_particles
            assert filled.contains_rigid is True

        filled.remove(filled.children[0])
        assert filled.contains_rigid is False
        assert filled.max_rigid_id is None

    def test_delete_body_semi_rigid(self, benzene):
        n_benzenes = 10
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4])
        filled.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')
        filled.remove(filled.children[0])

        assert filled.max_rigid_id == n_benzenes - 2
        assert len(list(filled.rigid_particles())) == (n_benzenes - 1) * 6

    def test_rigid_with_subcompounds1(self, rigid_benzene):
        compound = mb.Compound(subcompounds=rigid_benzene)

        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 0
        assert rigid_benzene.contains_rigid is True
        assert rigid_benzene.rigid_id is None
        assert rigid_benzene.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 12
        assert len(list(compound.rigid_particles(rigid_id=0))) == 12

    def test_rigid_with_subcompounds2(self, rigid_benzene):
        rigid_benzene2 = mb.clone(rigid_benzene)
        compound = mb.Compound(subcompounds=[rigid_benzene, rigid_benzene2])

        assert compound.max_rigid_id is 1
        assert rigid_benzene.max_rigid_id is 0
        assert rigid_benzene2.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 24
        assert len(list(compound.rigid_particles(rigid_id=0))) == 12
        assert len(list(compound.rigid_particles(rigid_id=1))) == 12

    def test_rigid_with_subcompounds3(self, benzene):
        benzene.label_rigid_bodies(rigid_particles='C')
        compound = mb.Compound(subcompounds=benzene)

        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 0
        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 6
        assert len(list(compound.rigid_particles(rigid_id=0))) == 6

    def test_rigid_with_subcompounds4(self, benzene):
        benzene.label_rigid_bodies(rigid_particles='C')
        benzene2 = mb.clone(benzene)
        compound = mb.Compound(subcompounds=[benzene, benzene2])

        assert compound.contains_rigid is True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 1
        assert benzene.contains_rigid is True
        assert benzene.rigid_id is None
        assert benzene.max_rigid_id is 0
        assert benzene2.contains_rigid is True
        assert benzene2.rigid_id is None
        assert benzene2.max_rigid_id is 1
        assert len(list(compound.rigid_particles())) == 12
        assert len(list(compound.rigid_particles(rigid_id=0))) == 6
        assert len(list(compound.rigid_particles(rigid_id=1))) == 6

    def test_rigid_with_subcompounds5(self, rigid_benzene):
        rigid_benzene2 = mb.clone(rigid_benzene)
        double = mb.Compound(subcompounds=[rigid_benzene, rigid_benzene2])
        double2 = mb.clone(double)
        compound = mb.Compound(subcompounds=[double, double2])

        assert compound.max_rigid_id is 3
        assert len(list(compound.rigid_particles())) == 48
        for rigid_id in range(4):
            assert len(list(compound.rigid_particles(rigid_id=rigid_id))) == 12

    def test_set_rigid_not_particle(self, benzene_from_parts):
        benzene_from_parts.label_rigid_bodies(rigid_particles=['C','H'])

        assert benzene_from_parts.contains_rigid is True
        assert benzene_from_parts.rigid_id is None
        assert benzene_from_parts.max_rigid_id is 0
        assert len(list(benzene_from_parts.rigid_particles())) == 12
        assert len(list(benzene_from_parts.rigid_particles(rigid_id=0))) == 12

    def test_manual_set_rigid_id(self, benzene):
        benzene[0].rigid_id = 0
        assert benzene.contains_rigid is True
        assert benzene[0].contains_rigid is False
        assert benzene.max_rigid_id is 0
        assert len(list(benzene.rigid_particles())) == 1

    def test_manual_set_rigid_id_error(self, benzene):
        with pytest.raises(AttributeError):
            benzene.rigid_id = 0

    def test_build_from_single_particle(self):
        compound = mb.Compound()
        compound.rigid_id = 0
        atom = mb.Compound(name='atom')
        atom.rigid_id = 0
        atom2 = mb.clone(atom)
        compound.add([atom, atom2], reset_rigid_ids=False)

        assert compound.contains_rigid == True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 2

    def test_build_from_single_particle2(self):
        compound = mb.Compound()
        compound.rigid_id = 0
        atom = mb.Compound(name='atom')
        atom.rigid_id = 0
        atom2 = mb.clone(atom)
        compound.add(atom)
        compound.add(atom2, reset_rigid_ids=False)

        assert compound.contains_rigid == True
        assert compound.rigid_id is None
        assert compound.max_rigid_id is 0
        assert len(list(compound.rigid_particles())) == 2
