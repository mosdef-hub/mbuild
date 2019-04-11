class Recipes(object):
    pass

recipes = Recipes()
from pkg_resources import iter_entry_points
available_methods = []
for entry_point in iter_entry_points(group='mbuild.plugins', name=None):
    setattr(recipes, entry_point.name, entry_point.load())
