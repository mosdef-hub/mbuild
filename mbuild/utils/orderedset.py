"""Ordered set module."""
from collections.abc import MutableSet


class OrderedSet(MutableSet):
    """An ordered set object with additional convenience methods.

    https://github.com/mosdef-hub/mbuild/issues/865

    Methods
    -------
    add
    discard
    index
    """

    def __init__(self, *args):
        self._data = {value: None for value in args}

    def __contains__(self, key):
        """Determine whether the set contains a key."""
        return key in self._data

    def __getitem__(self, value):
        """Get an item."""
        return list(self._data)[value]

    def __iter__(self):
        """Iterate through a copy."""
        return iter(list(self._data))

    def __len__(self):
        """Return the length."""
        return len(self._data)

    def add(self, value):
        """Add a value."""
        self._data[value] = None

    def discard(self, value):
        """Remove a value."""
        self._data.pop(value, None)

    def index(self, value):
        """Get the index of a value."""
        return list(self._data).index(value)
