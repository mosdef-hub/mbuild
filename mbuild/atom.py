from mbuild.mbase import MBase
from mbuild.part_mixin import PartMixin

__author__ = 'sallai'

from copy import deepcopy

import numpy as np

from mbuild.bond import Bond


class Atom(MBase, PartMixin):
    """Elementary container class - typically a leaf in the hierarchy.

    Note:
        -Also used as "ghost" particles in Ports.
        -Atoms can be added and substracted using +/- operators. The result is
        the addition or subtraction of the Atoms' cartesian coordinates.

    Attributes:
        kind (str): The kind of atom, usually the chemical element.
        pos (np.ndarray of floats): Cartesian coordinates of the atom.
        charge (float): Partial charge on the atom.
        parent (Compound): Compound to which the Atom belongs.
        referrers (set of Compounds): All Compounds that refer to this instance
            of Atom.
        bonds (set of Bond): Every Bond that the Atom is a part of.

    """
    __slots__ = ['kind', 'pos', 'charge', 'parent', 'referrers', 'bonds', 'uid']

    def __init__(self, kind, pos=None, charge=0.0):
        """Initialize an Atom.

        Args:
            kind (str): The kind of atom, usually the chemical element.
            pos (np.ndarray, optional): Cartesian coordinates of the atom.
            charge (float, optional): Partial charge on the atom.

        """
        super(Atom, self).__init__(kind, pos=pos, charge=charge)

        assert (isinstance(kind, basestring))

        if pos is None:
            pos = np.array([0, 0, 0])

        self.kind = kind
        self.pos = pos
        self.charge = charge
        # self.parent = None
        # self.referrers = set()
        self.bonds = set()


    def bonded_atoms(self, memo=dict()):
        """Return a list of atoms bonded to self. """
        for bond in self.bonds:
            bonded_atom = bond.other_atom(self)
            if id(bonded_atom) not in memo:
                memo[id(bonded_atom)] = bonded_atom
                bonded_atom.bonded_atoms(memo)
        return memo.values()

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

        # remember the topmost component being deepcopied
        if len(memo) == 0:
            memo[0] = self
        memo[id(self)] = newone

        # copy fields that don't need recursion
        newone.referrers = set()
        newone.bonds = set()
        newone.kind = deepcopy(self.kind, memo)
        newone.pos = deepcopy(self.pos, memo)
        newone.charge = deepcopy(self.charge, memo)

        # copy the parent of everybody, except the topmost compound being deepcopied
        if memo[0] == self or isinstance(memo[0], Bond):
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)

        return newone

