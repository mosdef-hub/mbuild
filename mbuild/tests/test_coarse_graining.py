import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn


class TestCoarseGraining(BaseTest):
    def test_hexane(self, hexane, propyl):
        particles = [propyl.__class__ ]
        cg = mb.coarse_grain(hexane, particle_classes=particles)
        assert cg.n_particles == 2
        assert cg.n_bonds == 1
        assert all(child.name.startswith(propyl.name)
                   for child in cg.children)

    def test_reverse_map_hexane(self, propane_aa):

        mapping_moieties = {'Propane': propane_aa}
        cg = mb.load(get_fn('hexane_cg.mol2')) 

        recovered = mb.reverse_map(cg, mapping_moieties, energy_minimize=False)
        assert recovered.n_particles == 20
        assert recovered.n_bonds == 19


