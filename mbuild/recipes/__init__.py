"""Entrypoints for mBuild recipe plugins."""

from importlib import metadata


class Recipes:
    """mBuild recipe object."""


recipes = Recipes()

entry_points = metadata.entry_points().select(group="mbuild.plugins")

available_methods = []
for entry_point in entry_points:
    setattr(recipes, entry_point.name, entry_point.load())
