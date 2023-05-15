"""Entrypoints for mBuild recipe plugins."""


class Recipes(object):
    """mBuild recipe object."""

    pass


recipes = Recipes()
from importlib import metadata

# from importlib import iter_entry_points
entry_points = metadata.entry_points()["mbuild.plugins"]
# entry_points = metadata.entry_points(group="mbuild.plugins", name=None)

available_methods = []
for entry_point in entry_points:
    setattr(recipes, entry_point.name, entry_point.load())
