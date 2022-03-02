import mbuild as mb
from mbuild.lib.moieties import CH3
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn


class TestMol2(BaseTest):
    def test_load_and_create(self):
        mb.load(get_fn("methyl.mol2"))

    def test_update_coordinates(self):
        methyl = CH3()
        methyl.update_coordinates(get_fn("methyl.mol2"))

    def test_save(self):
        methyl = mb.load(get_fn("methyl.mol2"))
        methyl.save(filename="methyl_out.mol2")

    def test_gmso_backend(self):
        pmd_silica_surface = mb.load(
            filename_or_object=get_fn("beta-cristobalite-expanded.mol2"),
            backend="parmed",
        )
        gmso_silica_surface = mb.load(
            filename_or_object=get_fn("beta-cristobalite-expanded.mol2"),
            backend="gmso",
        )

        assert pmd_silica_surface.n_particles == gmso_silica_surface.n_particles
        assert pmd_silica_surface.n_bonds == gmso_silica_surface.n_bonds

        element_set = set()
        for particle in gmso_silica_surface:
            element_set.add(particle.element)
        assert len(element_set) == 2
