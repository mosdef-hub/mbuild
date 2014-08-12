from mbuild.bond import Bond

__author__ = 'sallai'
from copy import deepcopy

import numpy as np

class Atom(object):
    """ """
    __slots__ = ['kind', 'pos', 'charge', 'parent', 'referrers', 'bonds', 'uid']

    def __init__(self, kind, pos=None, charge=0.0):
        assert (isinstance(kind, basestring))

        if pos is None:
            pos = np.array([0, 0, 0])

        self.kind = kind
        self.pos = pos
        self.charge = charge
        self.parent = None
        self.referrers = set()
        self.bonds = set()

    def ancestors(self):
        """
        Generate all ancestors of the Compound recursively
        :yield ancestors
        """
        yield self.parent
        if self.parent is not None:
            for a in self.parent.ancestors():
                yield a

    def __add__(self, other):
        if isinstance(other, Atom):
            other = other.pos
        return self.pos + other

    def __radd__(self, other):
        if isinstance(other, Atom):
            other = other.pos
        return self.pos + other

    def __sub__(self, other):
        if isinstance(other, Atom):
            other = other.pos
        return self.pos - other

    def __rsub__(self, other):
        if isinstance(other, Atom):
            other = other.pos
        return self.pos - other

    def __neg__(self):
        return -self.pos

    def __repr__(self):
        return "Atom{0}({1}, {2})".format(id(self), self.kind, self.pos)

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        if len(memo) == 0:
            memo[0] = self
        memo[id(self)] = newone

        newone.kind = deepcopy(self.kind, memo)
        newone.pos = deepcopy(self.pos, memo)
        newone.charge = deepcopy(self.charge, memo)

        # copy the parent of everybody, except the topmost compound being deepcopied
        if memo[0] == self or isinstance(memo[0], Bond):
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)

        newone.referrers = set()
        newone.bonds = set()

        return newone
