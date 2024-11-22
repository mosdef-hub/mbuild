"""Ordered set module."""

from collections.abc import MutableSet


class OrderedSet(MutableSet):
    """An ordered set object with additional convenience methods.

    Taken from code suggested by Vyas in
    https://github.com/mosdef-hub/mbuild/issues/865

    Methods
    -------
    add
    discard
    union
    intersection
    difference
    """

    def __init__(self, *args):
        self._data = {value: None for value in args}

    def __repr__(self):
        """Return the OrderedSet representation."""
        data = ", ".join(str(i) for i in self._data)
        return f"{self.__class__.__name__}({data})"

    def __contains__(self, key):
        """Determine whether the element `key` is in the set."""
        return key in self._data

    def __getitem__(self, value):
        """Get an item at index `value`."""
        return list(self._data)[value]

    def __iter__(self):
        """Iterate through the set."""
        return iter(self._data)

    def __len__(self):
        """Return the length."""
        return len(self._data)

    def add(self, value):
        """Add a value."""
        self._data[value] = None

    def discard(self, value):
        """Remove a value."""
        self._data.pop(value, None)

    def remove(self, value):
        """Remove a value. Alias for discard."""
        self.discard(value)

    def union(self, iterable):
        """Return the union of this set and an iterable."""
        new = OrderedSet()
        new._data = {
            **{i: None for i in self._data},
            **{i: None for i in iterable},
        }
        return new

    def intersection(self, iterable):
        """Return the intersection of this set and an iterable."""
        new = OrderedSet()
        new._data = {i: None for i in self._data if i in iterable}
        return new

    def difference(self, iterable):
        """Return the difference of this set and an iterable."""
        new = OrderedSet()
        new._data = {i: None for i in self._data if i not in iterable}
        return new
