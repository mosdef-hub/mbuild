from collections import OrderedDict
from copy import deepcopy

import numpy as np

from atom import Atom
from bond import Bond
from box import Box
from mbuild.periodic_kdtree import PeriodicCKDTree
from orderedset import OrderedSet


class Compound(object):
    """ """

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

        self.parent = None
        self.referrers = set()

    def add(self, new_obj, label=None, containment=True, replace=False,
            inherit_periodicity=True):
        """ """
        assert isinstance(new_obj, (Atom, Bond, Compound, list, tuple))
        if containment:
            # support batch add
            if isinstance(new_obj, (list, tuple)):
                for elem in new_obj:
                    assert(elem.parent is None)
                    self.add(elem)
                    elem.parent = self
                return
            assert(new_obj.parent is None)
            self.parts.add(new_obj)
            new_obj.parent = self

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
                new_obj.referrers.add(self)

    def remove(self, objs_to_remove):
        """
        Remove a part (atom, bond or component) from the compound by value
        :param obj: the part to remove
        """

        if not isinstance(objs_to_remove, (list, tuple, set)):
            objs_to_remove = set([objs_to_remove])

        if len(objs_to_remove) == 0:
            return

        intersection = objs_to_remove.intersection(self.parts)
        self.parts.difference_update(intersection)
        objs_to_remove.difference_update(intersection)

        for removed_part in intersection:
            removed_part.parent = None
            # remove labels in the hierarchy pointing to this part
            referrers_to_remove = set()
            for referrer in removed_part.referrers:
                if not removed_part in referrer.ancestors():
                    for label, referred_part in referrer.labels.items():
                        if referred_part is removed_part:
                            del referrer.labels[label]
                            # removed_part.referrers.remove(referrer)
                            referrers_to_remove.add(referrer)
            removed_part.referrers.difference_update(referrers_to_remove)

            # remove labels in this part pointing into the hierarchy
            labels_to_delete = []
            if isinstance(removed_part, Compound):
                for k, v in removed_part.labels.items():
                    if not removed_part in v.ancestors():
                        v.referrers.remove(removed_part)
                        # del removed_part.labels[k]
                        labels_to_delete.append(k)
            for k in labels_to_delete:
                del removed_part.labels[k]

            # if removing an atom, make sure to remove the bonds it's part of
            if isinstance(removed_part, Atom):
                for bond in removed_part.bonds:
                    if bond.parent is not None:
                        bond.parent.remove(bond)

        # remove it recursively from sub-components
        for part in self.parts:
            if isinstance(part, Compound) and len(objs_to_remove) > 0:
                part.remove(objs_to_remove)


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

    def referenced_ports(self):
        from mbuild.port import Port
        return [port for port in self.labels.values() if isinstance(port, Port)]

    def ancestors(self):
        """Generate all ancestors of the Compound recursively.

        Yields:
            ancestor (Compound): A Compound one or more levels higher in the
                hierarchy.

        """
        yield self.parent
        if self.parent is not None:
            for a in self.parent.ancestors():
                yield a

    def initAtomsByKind(self, kind='*'):
        # remember the hash of the parts dict at time of generating the atomListsByKind dict
        self.atomListsByKind_hash = self.parts.__hash__
        self.atomListsByKind = OrderedDict()

        # print "intializing atoms by kind dict"

        self.atomListsByKind['*'] = []

        for atom in self.atoms():

            self.atomListsByKind['*'].append(atom)

            if atom.kind not in self.atomListsByKind:
                self.atomListsByKind[atom.kind] = [atom]
            else:
                self.atomListsByKind[atom.kind].append(atom)

    def hasAtomListByKind(self, kind='*'):
        if not hasattr(self, 'atomListsByKind') or self.parts.__hash__ != self.atomListsByKind_hash:
            # print "nonexistent of outdated atomsListByKind"
            return False
        else:
            return True

    def getAtomListByKind(self, kind='*'):
        # use some precomputed data structures instead (memory vs. time tradeoff)
        if not self.hasAtomListByKind(kind):
            self.initAtomsByKind(kind)

        if kind in self.atomListsByKind:
            return self.atomListsByKind[kind]
        else:
            return []


    def initAtomKdTree(self, kind='*'):
            # check if atomKdTrees dict exists and is up-to-date
            if not hasattr(self, 'atomKdTrees') or self.parts.__hash__ != self.atomKdTrees_hash:
                # remember the hash of the bonds dict at time of generating the bondsByAtomKind dict
                self.atomKdTrees_hash = self.parts.__hash__
                self.atomKdTrees = dict()
                # print "intiializing atomKdTrees dict"

            #self.atomKdTrees[kind] = PeriodicCKDTree([atom.pos for atom in self.getAtomListByKind(kind)], bounds=self.periodicity)
            # host_atom_list = [atom for atom in self.atoms()]
            # host_atom_pos_list = [atom.pos for atom in host_atom_list]

            atom_list = self.getAtomListByKind(kind)
            atom_pos_list = [atom.pos for atom in atom_list]
            if len(atom_list) > 0:
                self.atomKdTrees[kind] = PeriodicCKDTree(atom_pos_list)
            else:
                self.atomKdTrees[kind] = None

    def hasAtomKdTree(self, kind='*'):
        if not hasattr(self, 'atomKdTrees') or self.parts.__hash__ != self.atomKdTrees_hash:
            return False

        if kind in self.atomKdTrees:
            return True

        return False

    def getAtomKdTree(self, kind='*'):
        return self.atomKdTrees[kind]

    def getAtomsInRange(self, point, radius, maxItems=10, kind='*'):

        # create kdtree if it's not yet there
        if not self.hasAtomKdTree(kind):
            self.initAtomKdTree(kind)

        if self.getAtomKdTree(kind) is None:
            return []

        distances, indices = self.getAtomKdTree(kind).query(point, maxItems)

        neighbors = []
        for index, distance in zip(indices, distances):
            if distance <= radius:
                al = self.getAtomListByKind(kind)
                neighbors.append(self.getAtomListByKind(kind)[index])
            else:
                break

        return neighbors


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

        return Box(mins=min_coords, maxes=max_coords)
        # return min_coords, max_coords, max_coords - min_coords

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        if len(memo) == 0:
            memo[0] = self
        memo[id(self)] = newone

        newone.kind = deepcopy(self.kind, memo)
        newone.periodicity = deepcopy(self.periodicity, memo)

        # copy the parent of everybody, except the topmost compound being deepcopied
        if memo[0] == self:
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)

        # copy parts, except bonds with atoms outside the hierarchy
        newone.parts = OrderedSet()
        for part in self.parts:
            if isinstance(part, Bond):
                if memo[0] in part.atom1.ancestors() and memo[0] in part.atom2.ancestors():
                    newone.parts.add(deepcopy(part,memo))
            else:
                newone.parts.add(deepcopy(part,memo))


        # copy labels, except bonds with atoms outside the hierarchy
        newone.labels = OrderedDict()
        for k, v in self.labels.items():
            if isinstance(v, Bond):
                if memo[0] in v.atom1.ancestors() and memo[0] in v.atom2.ancestors():
                    newone.labels[k] = deepcopy(v, memo)
                    newone.labels[k].referrers.add(newone)
            else:
                newone.labels[k] = deepcopy(v, memo)
                newone.labels[k].referrers.add(newone)


        # copy references that don't point out of the hierarchy
        newone.referrers = set()
        # for r in self.referrers:
        #     if memo[0] in r.ancestors():
        #         newone.referrers.add(r)

        return newone
