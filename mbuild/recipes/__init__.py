"""Entrypoints for mBuild recipe plugins."""

import sys
from importlib import metadata


class Recipes(object):
    """mBuild recipe object."""

    pass


recipes = Recipes()

if sys.version_info.minor >= 10:
    entry_points = metadata.entry_points().select(group="mbuild.plugins")
else:
    entry_points = metadata.entry_points()["mbuild.plugins"]

available_methods = []
for entry_point in entry_points:
    setattr(recipes, entry_point.name, entry_point.load())
