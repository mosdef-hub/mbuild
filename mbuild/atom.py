from copy import deepcopy

__author__ = 'sallai'

# import inspect
# import sys

from coordinate_transform import *


class Atom(object):
    __slots__ = ['kind', 'pos', 'charge', 'bonds', 'angles', 'dihedrals', 'parent']

    def __init__(self, kind, pos=None, charge=0):
        assert (isinstance(kind, basestring))

        if pos is None:
            pos = np.array([0, 0, 0])

        self.parent = None
        self.kind = kind
        self.pos = pos
        self.charge = charge
        self.bonds = set()
        self.angles = set()
        self.dihedrals = set()

    def distance(self, a2):
        return np.linalg.norm(self.pos - a2.pos)
        # return sqrt((self.pos[0] - a2.pos[0]) ** 2 + (self.pos[1] - a2.pos[1]) ** 2 + (self.pos[2] - a2.pos[2]) ** 2)

    def ancestors(self):
        ancestors = []

        current = self
        while current.parent is not None:
            ancestors.append(current.parent)
            current = current.parent

        ancestors.reverse()

        return ancestors

    def treeDistance(self, a2):

        a1 = self

        ancestors1 = a1.ancestors()
        ancestors2 = a2.ancestors()

        # print "ancestors1: " + str(ancestors1)
        # print "ancestors2: " + str(ancestors2)


        if ancestors1[0] != ancestors2[0]:
            # they are in separate component trees
            return float("inf")

        i = 0

        while len(ancestors1)>i and len(ancestors2)>i and ancestors1[i] == ancestors2[i]:
            i += 1

        # return len(ancestors1) + len(ancestors2) - 2*i

        treed = 0
        for c in ancestors1[i:]:
            treed += c.treeDistancePenalty

        for c in ancestors2[i:]:
            treed += c.treeDistancePenalty

        # print str(treed) + " " + str(len(ancestors1) + len(ancestors2) - 2*i)
        # assert(treed==len(ancestors1) + len(ancestors2) - 2*i)

        return treed

    def transform(self, T):
        self.pos = T.apply(self.pos)

    def __repr__(self):
        return "Atom" + str(id(self)) + "(" + self.kind + "," + str(self.pos) + ")"

    def __copy__(self):
        cls = self.__class__
        newone = cls.__new__(cls)

        newone.__dict__.update(self.__dict__)

        newone.parent = None

        return newone

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)

        memo[id(self)] = newone

        if id(self.parent) in memo:
            newone.parent = memo[id(self.parent)]
        else:
            newone.parent = None
        newone.kind = deepcopy(self.kind, memo)
        newone.pos = deepcopy(self.pos, memo)
        newone.charge = deepcopy(self.charge, memo)
        newone.bonds = deepcopy(self.bonds, memo)
        newone.angles = deepcopy(self.angles, memo)
        newone.dihedrals = deepcopy(self.dihedrals, memo)

        # for k, v in self.__dict__.items():
        #     setattr(newone, k, deepcopy(v, memo))

        return newone

