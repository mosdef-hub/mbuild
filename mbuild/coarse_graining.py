"""Tools for coarse-grained systems."""

from collections import OrderedDict
from copy import deepcopy

from mbuild.compound import Compound, clone
from mbuild.exceptions import MBuildError

__all__ = ["coarse_grain"]


def coarse_grain(real_thing, memo=None, particle_classes=None):
    """Coarse-graining function."""
    if memo is None:
        memo = OrderedDict()

    if particle_classes is None:
        particle_classes = []

    proxy = _create_proxy_compounds(real_thing, memo, particle_classes)
    _create_proxy_bonds(real_thing, memo, particle_classes)
    _create_proxy_labels(real_thing, memo)

    return proxy


class Proxy(Compound):
    """Proxy class."""

    def __init__(self, compound):
        if compound.name == "G":
            name = "G"
        else:
            name = compound.name + "_PROXY"
        super(Proxy, self).__init__(name=name)

        self.wrapped = compound
        self.children = None
        self.labels = None
        self.parent = None
        self.referrers = set()
        self.index = None
        self.bond_graph = None

    def proxy_for(self):
        """Get the proxy."""
        if hasattr(self.wrapped, "wrapped"):
            return self.wrapped.proxy_for()
        else:
            return self.wrapped.__class__

    @property
    def pos(self):
        """Get the wrapped position."""
        return self.wrapped.pos

    @pos.setter
    def pos(self, value):
        self.wrapped.pos = value

    def __getattr__(self, attr):
        """Get the attribute from the class."""
        return getattr(self.wrapped, attr)

    def _clone(self, clone_of=None, root_container=None):
        """Clone a Proxy.

        A faster alternative to deepcopying. Does not resolve circular
        dependencies. This should be safe provided you never try to add the
        top of a Compound hierarchy to a sub-Compound. Clones compound
        hierarchy only, not the bonds.
        """
        if root_container is None:
            root_container = self
        if clone_of is None:
            clone_of = dict()

        # If this compound has already been cloned, return that.
        if self in clone_of:
            return clone_of[self]

        # Otherwise we make a new clone.
        cls = self.__class__
        newone = cls.__new__(cls)

        # Remember that we're cloning the new one of self.
        clone_of[self] = newone

        newone.name = deepcopy(self.name)
        newone.wrapped = clone(self.wrapped)

        if hasattr(self, "index"):
            newone.index = deepcopy(self.index)

        if self.children is None:
            newone.children = None
        else:
            newone.children = list()
        # Parent should be None initially.
        newone.parent = None
        newone.labels = OrderedDict()
        newone.referrers = set()
        newone.bond_graph = None

        # Add children to clone.
        if self.children:
            for child in self.children:
                newchild = child._clone(clone_of, root_container)
                newone.children.append(newchild)
                newchild.parent = newone

        # Copy labels, except bonds with atoms outside the hierarchy.
        if self.labels:
            for label, compound in self.labels.items():
                if not isinstance(compound, list):
                    newone.labels[label] = compound._clone(
                        clone_of, root_container
                    )
                    compound.referrers.add(clone_of[compound])
                else:
                    # compound is a list of compounds, so we create an empty
                    # list, and add the clones of the original list elements.
                    newone.labels[label] = []
                    for subpart in compound:
                        newone.labels[label].append(
                            subpart._clone(clone_of, root_container)
                        )
                        # Referrers must have been handled already, or the will
                        # be handled

        return newone

    def _clone_bonds(self, clone_of=None):
        newone = clone_of[self]
        for c1, c2 in self.bonds():
            try:
                newone.add_bond((clone_of[c1], clone_of[c2]))
            except KeyError:
                raise MBuildError(
                    "Cloning failed. Compound contains bonds to "
                    "Particles outside of its containment hierarchy."
                )


def is_leaf(what):
    return hasattr(what, "parts") and not what.children


def _create_proxy_compounds(real_thing, memo, particle_classes):
    proxy = Proxy(real_thing)
    memo[real_thing] = proxy

    if not type(real_thing) in particle_classes:
        if not is_leaf(real_thing):  # Recurse only if it has parts.
            # Recursively create proxies for parts.
            for part in real_thing.children:
                part_proxy = _create_proxy_compounds(
                    part, memo, particle_classes
                )
                proxy.add(part_proxy)

    return proxy


def _proxy_of(real_thing, memo):
    if real_thing in memo:
        return memo[real_thing]
    else:
        return _proxy_of(real_thing.parent, memo)


def _create_proxy_bonds(real_thing, memo, leaf_classes):
    proxy = memo[real_thing]

    if type(real_thing) in leaf_classes or is_leaf(real_thing):
        # It is a leaf of the proxy, so we don't recurse.
        pass
    else:
        for part in real_thing.children:
            _create_proxy_bonds(part, memo, leaf_classes)

    # Check if there's a contained bond that needs to be added to the proxy.
    if hasattr(real_thing, "bonds"):
        for a1, a2 in real_thing.bonds():
            pa1 = _proxy_of(a1, memo)
            pa2 = _proxy_of(a2, memo)
            if pa1 != pa2:  # Do not add internal bonds.
                proxy.add_bond((pa1, pa2))


def _create_proxy_labels(real_thing, memo):
    if not is_leaf(real_thing):
        for label, part in real_thing.labels.items():
            if isinstance(part, list):
                continue
            if part in memo:
                memo[real_thing].labels[label] = memo[part]

        for part in real_thing.children:
            _create_proxy_labels(part, memo)
