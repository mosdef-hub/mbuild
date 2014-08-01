from collections import OrderedDict, defaultdict
from copy import deepcopy
from orderedset import OrderedSet
from warnings import warn

import numpy as np
from numpy.linalg import norm

from atom import Atom
from bond import Bond


class Compound(object):
    """ """
    __slots__ = ['kind', 'periodicity', 'parts', 'labels']

    def __init__(self, kind=None, periodicity=None):
        # set the kind (defaults to classname)
        if kind:
            self.kind = kind
        else:
            self.kind = self.__class__.__name__

        if not periodicity:
            periodicity = np.array([0.0, 0.0, 0.0])
        self.periodicity = periodicity

        # contains children (compounds or atoms)
        self.parts = OrderedSet()
        # label to compound/atom mappings -- need not be in parts
        self.labels = OrderedDict()

    def add(self, new_obj, label=None, containment=True, replace=False,
            inherit_periodicity=True):
        """ """
        assert isinstance(new_obj, (Atom, Bond, Compound, list, tuple))
        if containment:
            # support batch add
            if isinstance(new_obj, (list, tuple)):
                for elem in new_obj:
                    self.add(elem)
                return
            self.parts.add(new_obj)

        if (inherit_periodicity
            and isinstance(new_obj, Compound)
            and np.any(new_obj.periodicity)):
            self.periodicity = new_obj.periodicity

        # add new_obj to labels
        assert isinstance(new_obj, (Atom, Bond, Compound))

        if not containment and label is None:
            label = '_'+new_obj.__class__.__name__+'_'+str(id(new_obj))

        if label is not None:
            if not replace and label in self.labels:
                raise Exception("Label {0} already exists in {1}".format(label, self))
            else:
                self.labels[label] = new_obj


    def remove(self, obj):
        """
        Remove a part (atom, bond or component) from the compound by value
        :param obj: the part to remove
        """

        # remove reference to obj
        for k in self.labels:
            if self.labels[k] == obj:
                del self.labels[k]

        # remove containment
        if obj in self.parts:
            self.parts.discard(obj)

        # remove it recursively from subcomponents
        for part in self.parts:
            if isinstance(part, Compound):
                part.remove(obj)

    def __getattr__(self, attr):
        if attr in self.labels:
            return self.labels[attr]
        else:
            raise AttributeError

    def _yield_parts(self, part_type):
        """
        Generate all atoms of the Compound recursively
        :return: label - atom pairs
        """
        for part in self.parts:
            # add local atoms
            if isinstance(part, part_type):
                yield part
                # add atoms in sub-components recursively
            if isinstance(part, Compound):
                for subpart in part._yield_parts(part_type):
                    yield subpart

    def _find_parent_of(self, obj):
        # check if self is the parent
        if obj in self.parts:
            return self

        # recursively search for parent in parts
        for part in self.parts:
            if isinstance(part, Compound):
                parent = part._find_parent_of(obj)
                if not parent is None:
                    return parent

        # no parent found
        return None

    def atoms(self):
        return self._yield_parts(Atom)

    def bonds(self):
        return self._yield_parts(Bond)

    def min_periodic_distance(self, x0, x1):
        """Vectorized distance calculation considering minimum image
        """
        d = np.abs(x0 - x1)
        d = np.where(d > 0.5 * self.periodicity, self.periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

    def boundingbox(self, excludeG=True):
        """
        Compute the bounding box of the compound
        :rtype : (minx, miny, minz), (maxx, maxy, maxz)
        """
        minx = float('inf')
        miny = float('inf')
        minz = float('inf')
        maxx = float('-inf')
        maxy = float('-inf')
        maxz = float('-inf')

        for a in self.atoms():
            if excludeG and a.kind == 'G':
                continue
            if a.pos[0] < minx:
                minx = a.pos[0]
            if a.pos[0] > maxx:
                maxx = a.pos[0]
            if a.pos[1] < miny:
                miny = a.pos[1]
            if a.pos[1] > maxy:
                maxy = a.pos[1]
            if a.pos[2] < minz:
                minz = a.pos[2]
            if a.pos[2] > maxz:
                maxz = a.pos[2]

        min_coords = np.array([minx, miny, minz])
        max_coords = np.array([maxx, maxy, maxz])

        return min_coords, max_coords, max_coords - min_coords

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        memo[id(self)] = newone

        newone.kind = deepcopy(self.kind, memo)
        newone.periodicity = deepcopy(self.periodicity, memo)
        newone.parts = deepcopy(self.parts, memo)
        newone.labels = deepcopy(self.labels, memo)

        return newone


