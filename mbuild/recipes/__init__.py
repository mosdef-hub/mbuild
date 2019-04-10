class Recipes(object):
    pass

Recipes = Recipes()
from pkg_resources import iter_entry_points
available_methods = []
for entry_point in iter_entry_points(group='mbuild.plugins', name=None):
    setattr(Recipes, entry_point.name, entry_point.load())
