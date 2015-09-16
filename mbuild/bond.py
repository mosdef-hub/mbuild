# from copy import deepcopy
#
# import numpy as np
#
# from mbuild.part import Part
#
#
# class Bond(Part):
#     """Connection between two Atoms.
#
#     Attributes
#     ----------
#     atom1 : mb.Atom
#         First Atom in the bond.
#     atom2 : mb.Atom
#         Second Atom in the bond.
#     parent : mb.Compound
#         Compound to which the Bond belongs.
#     """
#     __slots__ = ['_atom1', '_atom2', 'kind', 'parent', 'referrers']
#
#     def __init__(self, atom1, atom2, kind=None):
#         super(Bond, self).__init__()
#         assert(not atom1 == atom2)
#
#         # If a Port is used to initialize a Bond, the Atom that the Port is
#         # anchored to will be used to create the Bond.
#         from mbuild.port import Port
#         if isinstance(atom1, Port):
#             atom1 = atom1.anchor
#         if isinstance(atom2, Port):
#             atom2 = atom2.anchor
#         self._atom1 = atom1
#         self._atom2 = atom2
#
#         if kind is not None:
#             self.kind = kind
#         else:
#             self.kind = '{0}-{1}'.format(atom1.name, atom2.name)
#
#         # Ensure Atoms in Bond know about the Bond.
#         if atom1.attached_bonds is None:
#             atom1.attached_bonds = set()
#         atom1.attached_bonds.add(self)
#
#         if atom2.attached_bonds is None:
#             atom2.attached_bonds = set()
#         atom2.attached_bonds.add(self)
#
#     @property
#     def atom1(self):
#         return self._atom1
#
#     @atom1.setter
#     def atom1(self, atom):
#         self._atom1 = atom
#
#     @property
#     def atom2(self):
#         return self._atom2
#
#     @atom2.setter
#     def atom2(self, atom):
#         self._atom2 = atom
#
#     def other_atom(self, atom):
#         """Returns the other Atom in the Bond. """
#         if self.atom1 is atom:
#             return self.atom2
#         elif self.atom2 is atom:
#             return self.atom1
#
#     def length(self, periodicity=None):
#         """Calculate the bond length considering minimum image. """
#         dist = np.abs(self.atom1.pos - self.atom2.pos)
#         if periodicity is None:
#             return np.sqrt((dist ** 2).sum(axis=-1))
#         else:
#             dist = np.where(dist > 0.5 * periodicity, periodicity - dist, dist)
#             return np.sqrt((dist ** 2).sum(axis=-1))
#
#     def has_both_atoms_within(self, root_compound):
#         """Both atoms are within the containment hierarchy of a Compound. """
#         return (root_compound in self.atom1.ancestors() or
#                 root_compound in self.atom2.ancestors())
#
#     def has_atoms_outside_of(self, root_compound):
#         """Has atoms outside of the containment hierarchy of a Compound."""
#         return (root_compound not in self.atom1.ancestors() or
#                 root_compound not in self.atom2.ancestors())
#
#     def __repr__(self):
#         return "Bond{0}({1}, {2})".format(id(self), self.atom1, self.atom2)
#
#     def _clone(self, clone_of=None, root_container=None):
#         from mbuild import clone
#
#         if clone_of is None:
#             clone_of = dict()
#
#         # If this bond has already been cloned, return that.
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
#         newone.kind = self.kind
#         newone.referrers = set()
#
#         # Do the rest recursively.
#         newone.atom1 = clone(self.atom1, clone_of, root_container)
#         newone.atom2 = clone(self.atom2, clone_of, root_container)
#         newone.atom1.attached_bonds.add(newone)
#         newone.atom2.attached_bonds.add(newone)
#
#         # We set newone.parent in compound.
#         return newone
#
#     def __deepcopy__(self, memo):
#         cls = self.__class__
#         newone = cls.__new__(cls)
#
#         # Remember the topmost component being deepcopied.
#         if len(memo) == 0:
#             print('bond is root of deepcopy')
#             memo[0] = self
#         memo[id(self)] = newone
#
#         # Copy fields that don't need recursion.
#         newone.kind = self.kind
#         newone.referrers = set()
#
#         # Do the rest recursively.
#         newone.atom1 = deepcopy(self.atom1, memo)
#         newone.atom2 = deepcopy(self.atom2, memo)
#         newone.atom1.attached_bonds.add(newone)
#         newone.atom2.attached_bonds.add(newone)
#
#         # Copy parents, except the topmost compound being deepcopied.
#         if memo[0] == self:
#             newone.parent = None
#         else:
#             newone.parent = deepcopy(self.parent, memo)
#
#         return newone
