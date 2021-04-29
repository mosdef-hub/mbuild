import mbuild as mb
from mbuild.lib.moieties import CH3
from mbuild.utils.io import get_fn
from mbuild.tests.base_test import BaseTest


class TestMol2(BaseTest):

    def test_load_and_create(self):
        mb.load(get_fn('methyl.mol2'))

    def test_update_coordinates(self):
        methyl = CH3()
        methyl.update_coordinates(get_fn('methyl.mol2'))

    def test_save(self):
        methyl = mb.load(get_fn('methyl.mol2'))
        methyl.save(filename='methyl_out.mol2')

