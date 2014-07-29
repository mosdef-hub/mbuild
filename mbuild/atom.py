__author__ = 'sallai'
from copy import deepcopy

import numpy as np

class Atom(object):
    """ """
    __slots__ = ['kind', 'pos', 'charge']

    def __init__(self, kind, pos=None, charge=0.0):
        assert (isinstance(kind, basestring))

        if pos is None:
            pos = np.array([0, 0, 0])

        self.kind = kind
        self.pos = pos
        self.charge = charge

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
        memo[id(self)] = newone

        newone.kind = deepcopy(self.kind, memo)
        newone.pos = deepcopy(self.pos, memo)
        newone.charge = deepcopy(self.charge, memo)
        return newone
