import pytest

from mbuild.tests.base_test import BaseTest


class TestPlugins(BaseTest):
    def test_basic_import(self):
        import mbuild.recipes

    @pytest.mark.parametrize(
        'recipe_name',
        ['Monolayer', 'Polymer', 'SilicaInterface', 'TiledCompound'],
    )
    def test_recipes_contents(self, recipe_name):
        import mbuild.recipes
        recipe_name in dir(mbuild.recipes)
