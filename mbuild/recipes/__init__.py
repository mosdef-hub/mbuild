"""Entrypoints for mBuild recipe plugins."""


class Recipes(object):
    """mBuild recipe object."""

    pass


recipes = Recipes()
from importlib import metadata

entry_points = metadata.entry_points()["mbuild.plugins"]

available_methods = []
for entry_point in entry_points:
    setattr(recipes, entry_point.name, entry_point.load())
