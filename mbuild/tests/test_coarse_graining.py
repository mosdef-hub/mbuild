import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestCoarseGraining(BaseTest):
    def test_hexane(self, hexane, propyl):
        particles = [propyl.__class__]
        cg = mb.coarse_grain(hexane, particle_classes=particles)
        assert cg.n_particles == 2
        assert cg.n_bonds == 1
        assert all(child.name.startswith(propyl.name) for child in cg.children)

        cg_clone = mb.clone(cg)
        assert cg_clone.n_particles == 2
        assert cg_clone.n_bonds == 1
        assert all(child.name.startswith(propyl.name) for child in cg_clone.children)
        assert cg_clone.wrapped.n_particles == 20
        assert cg_clone.wrapped.n_bonds == 19
