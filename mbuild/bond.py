
import numpy as np

from mbuild.part import Part


class Bond(Part):
    """Connection between two Atoms.

    Attributes
    ----------
    atom1 : mb.Atom
        First Atom in the bond.
    atom2 : mb.Atom
        Second Atom in the bond.
    parent : mb.Compound
        Compound to which the Bond belongs.
    """
    __slots__ = ['_atom1', '_atom2', 'kind', 'parent', 'referrers']

    def __init__(self, atom1, atom2, kind=None):
        super(Bond, self).__init__()
        assert(not atom1 == atom2)

        # If a Port is used to initialize a Bond, the Atom that the Port is
        # anchored to will be used to create the Bond.
        from mbuild.port import Port
        if isinstance(atom1, Port):
            atom1 = atom1.anchor
        if isinstance(atom2, Port):
            atom2 = atom2.anchor
        self._atom1 = atom1
        self._atom2 = atom2

        if kind is not None:
            self.kind = kind
        else:
            self.kind = '{0}-{1}'.format(atom1.name, atom2.name)

        # Ensure Atoms in Bond know about the Bond.
        atom1.bonds.add(self)
        atom2.bonds.add(self)

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    def other_atom(self, atom):
        """Returns the other Atom in the Bond. """
        if self._atom1 is atom:
            return self._atom2
        elif self._atom2 is atom:
            return self._atom1

    def length(self, periodicity=None):
        """Calculate the bond length considering minimum image. """
        d = np.abs(self.atom1 - self.atom2)
        if periodicity is None:
            return np.sqrt((d ** 2).sum(axis=-1))
        else:
            d = np.where(d > 0.5 * periodicity, periodicity - d, d)
            return np.sqrt((d ** 2).sum(axis=-1))

    def has_both_atoms_within(self, root_compound):
        """Check both atoms are within the containment hierarchy of a root_compound."""
        return root_compound in self.atom1.ancestors() or root_compound in self.atom2.ancestors()

    def has_atoms_outside_of(self, root_compound):
        """Check either (or both) atoms are outside of the containment hierarchy of a root_compound."""
        return root_compound not in self.atom1.ancestors() or root_compound not in self.atom2.ancestors()

    def __repr__(self):
        return "Bond{0}({1}, {2})".format(id(self), self.atom1, self.atom2)

    def _clone(self, clone_of=None, root_container=None):
        from mbuild import clone
        # create the clone_of dict if it's None
        if not clone_of:
            clone_of=dict()

        # if this bond has been cloned, return it
        if self in clone_of:
            return clone_of[self]

        # else we make a new clone

        cls = self.__class__
        newone = cls.__new__(cls)

        # remember that we're cloning the new one of of self
        clone_of[self] = newone

        # Copy fields that don't need recursion.
        newone.kind = self.kind
        newone.referrers = set()

        # Do the rest recursively.
        newone._atom1 = clone(self.atom1, clone_of, root_container)
        newone._atom2 = clone(self.atom2, clone_of, root_container)
        newone._atom1.bonds.add(newone)
        newone._atom2.bonds.add(newone)

        # we set newone.parent in compound

        return newone


    def __deepcopy__(self, memo):
        from copy import deepcopy

        cls = self.__class__
        newone = cls.__new__(cls)

        # Remember the topmost component being deepcopied.
        if len(memo) == 0:
            print('bond is root of deepcopy')
            memo[0] = self
        memo[id(self)] = newone

        # Copy fields that don't need recursion.
        newone.kind = self.kind
        newone.referrers = set()

        # Do the rest recursively.
        newone._atom1 = deepcopy(self.atom1, memo)
        newone._atom2 = deepcopy(self.atom2, memo)
        newone._atom1.bonds.add(newone)
        newone._atom2.bonds.add(newone)

        # Copy parents, except the topmost compound being deepcopied.
        if memo[0] == self:
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)

        return newone
