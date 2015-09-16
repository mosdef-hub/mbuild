from mbuild.compound import Atom
__all__ = ['Atom']

# import numpy as np
#
# from mbuild.bond import Bond
# from mbuild.part import Part


# class Atom(Part):
#     """Elementary container class - typically a leaf in the hierarchy.
#
#     Notes
#     -----
#     Atoms are also used as "ghost" particles in Ports; such atoms are named 'G'.
#
#     Atoms can be added and substracted using +/- operators. The result is
#     the addition or subtraction of the Atoms' cartesian coordinates.
#
#     Attributes
#     ----------
#     name : str
#         The name of the atom, usually the chemical element.
#     pos : np.ndarray, shape=(3,), dtype=float
#         Cartesian coordinates of the atom.
#     charge : float
#         Partial charge on the atom.
#     parent : mb.Compound
#         Compound to which the Atom belongs.
#     referrers : set of mb.Compounds
#         All Compounds that refer to this instance of Atom.
#     bonds : set of mb.Bonds
#         Every Bond that the Atom is a part of.
#
#     """
#     # __slots__ = ['index', 'name', 'pos', 'charge', 'parent', 'referrers', 'bonds']
#
#     def __init__(self, name, pos=None, charge=0.0):
#         super(Atom, self).__init__()
#
#         if pos is None:
#             pos = np.array([0, 0, 0], dtype=float)
#
#         self.index = None  # Only used for specific purposes, e.g. TiledCompound
#         self.name = name
#         self.pos = np.asarray(pos, dtype=float)
#         self.charge = charge
#         self.bonds = set()
#
#     @property
#     def xyz_with_ports(self):
#         """Work around to make coordinate transforms work on Atoms. """
#         return self.pos
#
#     @property
#     def center(self):
#         """Work around to make coordinate transforms work on Atoms. """
#         return self.pos
#
#     @property
#     def neighbors(self):
#         """Return a list of all neighboring Atoms. """
#         return [bond.other_atom(self) for bond in self.bonds]
#
#     @property
#     def n_bonds(self):
#         return len(self.bonds)
#
#     def __add__(self, other):
#         if isinstance(other, Atom):
#             other = other.pos
#         return self.pos + other
#
#     def __radd__(self, other):
#         if isinstance(other, Atom):
#             other = other.pos
#         return self.pos + other
#
#     def __sub__(self, other):
#         if isinstance(other, Atom):
#             other = other.pos
#         return self.pos - other
#
#     def __rsub__(self, other):
#         if isinstance(other, Atom):
#             other = other.pos
#         return self.pos - other
#
#     def __neg__(self):
#         return -self.pos
#
#     def __repr__(self):
#         return "Atom{0}({1}, {2})".format(id(self), self.name, self.pos)
#
#     def _clone(self, clone_of=None, root_container=None):
#         if clone_of is None:
#             clone_of = dict()
#
#         # If this atom has already been cloned, return that.
#         if self in clone_of:
#             return clone_of[self]
#
#         # Otherwise, we make a new clone.
#         cls = self.__class__
#         newone = cls.__new__(cls)
#
#         # Remember that we're cloning the new one of self.
#         clone_of[self] = newone
#
#         # Copy fields that don't need recursion.
#         newone.referrers = set()  # We fill up the referrers set in compound.
#         newone.bonds = set()  # We fill up the bonds set in bonds.
#
#         # Do the rest recursively.
#         newone.index = self.index
#         newone.name = self.name
#         newone.pos = self.pos
#         newone.charge = self.charge
#
#         # We set newone.parent in compound.
#         return newone
#
#     def __deepcopy__(self, memo):
#         from copy import deepcopy
#
#         cls = self.__class__
#         newone = cls.__new__(cls)
#
#         # Remember the topmost component being deepcopied.
#         if len(memo) == 0:
#             memo[0] = self
#         memo[id(self)] = newone
#
#         # Copy fields that don't need recursion.
#         newone.referrers = set()
#         newone.bonds = set()
#
#         # Do the rest recursively.
#         newone.index = deepcopy(self.index, memo)
#         newone.name = deepcopy(self.name, memo)
#         newone.pos = deepcopy(self.pos, memo)
#         newone.charge = deepcopy(self.charge, memo)
#
#         # Copy parents, except the topmost compound being deepcopied.
#         if memo[0] == self or isinstance(memo[0], Bond):
#             newone.parent = None
#         else:
#             newone.parent = deepcopy(self.parent, memo)
#
#         return newone
