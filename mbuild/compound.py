from collections import OrderedDict
from copy import deepcopy

import numpy as np

from atom import Atom
from bond import Bond
from box import Box
from periodic_kdtree import PeriodicCKDTree
from orderedset import OrderedSet


class Compound(object):
    """A building block in the mBuild hierarchy.
   
    *** add eloquent description here ***

    Attributes:
        kind (str): The type of Compound.
        periodicity (np.ndarray): The periodic distance of the
                Compound in the x, y and z directions.
        parts (OrderedSet): All Atom, Bond and Compound instances that this
            Compound is a parent of.
        labels (OrderedDict): Strings referring to parts. This is primarily a
            convenience functionality.
        parent (Compound): Compound to which this Compound belongs.
        referrers (set of Compound): Compounds with labels referring to this
            Compound.

    """

    def __init__(self, kind=None, periodicity=None):
        """Construct a new Compound. 
        
        Args:
            kind (str, optional): The type of Compound.
            periodicity (np.ndarray, optional): The periodic distance of the
                Compound in the x, y and z directions.

        """
        # Set kind to classname if not specified.
        if kind:
            self.kind = kind
        else:
            self.kind = self.__class__.__name__

        # A periodocity of zero in any direction is treated as non-periodic.
        if not periodicity:
            periodicity = np.array([0.0, 0.0, 0.0])
        self.periodicity = periodicity

        # Contains all child parts. Parts can be Atom, Bond or Compound.
        self.parts = OrderedSet()

        # Labels to Compound/Atom mappings. These do not necessarily need not
        # be in self.parts.
        self.labels = OrderedDict()

        # The parent compound that contains this compound (can be None if this compound is the root of the
        # containment hierarchy.)
        self.parent = None

        # The set of other compounds that reference this compound with labels
        self.referrers = set()

        self._label_counter = dict()


    def atoms(self):
        return self._yield_parts(Atom)

    def bonds(self):
        return self._yield_parts(Bond)

    def _yield_parts(self, part_type):
        """Generate all parts of a specified type in the Compound recursively.

        Args:
            part_type (Atom, Bond, Compound): The type of parts to yield.

        Yields:
            part/subpart (Atom, Bond, Compound): A part in the hierarchy
                matching the specified type.

        """
        for part in self.parts:
            # Parts local to the current Compound.
            if isinstance(part, part_type):
                yield part
            # Parts in Compounds further down the hierarchy.
            if isinstance(part, Compound):
                for subpart in part._yield_parts(part_type):
                    yield subpart

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
            for ancestor in self.parent.ancestors():
                yield ancestor

    def add(self, new_part, label=None, containment=True, replace=False,
            inherit_periodicity=True):
        """Add a part to the Compound.

        Note: 
            This does not necessarily add the part to self.parts but may
            instead be used to add a reference to the part to self.labels. See
            'containment' argument.
       
        Args:
            new_part (Atom, Bond or Compound): The object to be added to this
                Compound.
            label (str, optional): A descriptive string for the part.
            containment (bool, optional):
            replace (bool, optional):
            inherit_periodicity (bool, optional):

        """
        assert isinstance(new_part, (Atom, Bond, Compound, list, tuple, set))
        if containment:
            # Support batch add via lists, tuples and sets.
            if isinstance(new_part, (list, tuple, set)):
                for elem in new_part:
                    assert(elem.parent is None)
                    self.add(elem)
                    elem.parent = self
                return
            # Is the following assertion necessary? Seems to be handled above.
            assert(new_part.parent is None)  
            self.parts.add(new_part)
            new_part.parent = self
        
        if (inherit_periodicity
            and isinstance(new_part, Compound)
            and np.any(new_part.periodicity)):
            self.periodicity = new_part.periodicity

        # Add new_part to labels. Does not currently support batch add.
        assert isinstance(new_part, (Atom, Bond, Compound))
        if not containment and label is None:
            label = '_{0}_{1}'.format(new_part.__class__.__name__, id(new_part))

        if label is not None:
            if label.endswith("[$]"):
                label = label[:-3]
                if not label in self.labels:
                    self.labels[label] = []

                label_pattern = label + "[{}]"

                count = len(self.labels[label])
                self.labels[label].append(new_part)
                label = label_pattern.format(count)

            if not replace and label in self.labels:
                raise Exception("Label {0} already exists in {1}".format(label, self))
            else:
                self.labels[label] = new_part
                new_part.referrers.add(self)

    def remove(self, objs_to_remove):
        """Remove a part (Atom, Bond or Compound) from the Compound by value.

        Args:
            objs_to_remove (set of parts): All objects to be removed from the
                hierarchy. If this is not a set, it will be cast to one to
                remove duplicates.
        """
        if not isinstance(objs_to_remove, (list, tuple, set)):
            objs_to_remove = [objs_to_remove]

        objs_to_remove = set(objs_to_remove)

        if len(objs_to_remove) == 0:
            # when does this occurr? 
            return

        intersection = objs_to_remove.intersection(self.parts)
        self.parts.difference_update(intersection)
        objs_to_remove.difference_update(intersection)

        for removed_part in intersection:
            removed_part.parent = None
            # Remove labels in the hierarchy pointing to this part.
            referrers_to_remove = set()
            for referrer in removed_part.referrers:
                if not removed_part in referrer.ancestors():
                    for label, referred_part in referrer.labels.items():
                        if referred_part is removed_part:
                            del referrer.labels[label]
                            referrers_to_remove.add(referrer)
            removed_part.referrers.difference_update(referrers_to_remove)

            # Remove labels in this part pointing into the hierarchy.
            labels_to_delete = []
            if isinstance(removed_part, Compound):
                for label, part in removed_part.labels.items():
                    if not removed_part in part.ancestors():
                        part.referrers.remove(removed_part)
                        labels_to_delete.append(label)
            for label in labels_to_delete:
                del removed_part.labels[label]

            # If removing an atom, make sure to remove the bonds it's part of.
            if isinstance(removed_part, Atom):
                for bond in removed_part.bonds:
                    if bond.parent is not None:
                        bond.parent.remove(bond)

        # Remove the part recursively from sub-components.
        for part in self.parts:
            if isinstance(part, Compound) and len(objs_to_remove) > 0:
                part.remove(objs_to_remove)

        # TODO: remove lists from labels if the list becomes empty, remove items from lists, as well

    def __getattr__(self, attr):
        if attr in self.labels:
            return self.labels[attr]
        else:
            raise AttributeError("'{}' object has no attribute '{}'".format(self.__class__.__name__, attr))

    def atom_list_by_kind(self, kind='*', excludeG=False, with_id_to_idx_mapping=False):
        list = []
        id_to_idx = dict()
        idx = 0
        for atom in self.atoms():
            if not (excludeG and atom.kind == "G"):
                if kind == '*':
                    list.append(atom)
                    id_to_idx[id(atom)] = idx
                    idx += 1
                elif atom.kind == kind:
                    list.append(atom)
                    id_to_idx[id(atom)] = idx
                    idx += 1

        if with_id_to_idx_mapping:
            return list, id_to_idx
        else:
            return list

    def bond_list_by_kind(self, kind='*', with_id_to_idx_mapping=False):
        list = []
        id_to_idx = dict()

        idx = 0
        for bond in self.bonds():
            if kind == '*':
                list.append(bond)
                id_to_idx[id(bond)] = idx
                idx += 1

            elif bond.kind == kind:
                list.append(bond)
                id_to_idx[id(bond)] = idx
                idx += 1

        if with_id_to_idx_mapping:
            return list, id_to_idx
        else:
            return list


    def min_periodic_distance(self, x0, x1):
        """Vectorized distance calculation considering minimum image. """
        d = np.abs(x0 - x1)
        d = np.where(d > 0.5 * self.periodicity, self.periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

    # def _init_atom_kdtree(self):
    #         atom_list = self.atom_list_by_kind('*', excludeG=True)
    #         atom_pos_list = [atom.pos for atom in atom_list]
    #         if len(atom_list) > 0:
    #             self._atom_kdtree = PeriodicCKDTree(atom_pos_list)
    #         else:
    #             self._atom_kdtree = None
    #
    # def atoms_in_range(self, point, radius, maxItems=10):
    #
    #     # create kdtree if it's not yet there
    #     if not hasattr(self,'_atom_kdtree'):
    #         self._init_atom_kdtree()
    #
    #     if self._atom_kdtree is None:
    #         return []
    #
    #     distances, indices = self._atom_kdtree.query(point, maxItems)
    #
    #     neighbors = []
    #     for index, distance in zip(indices, distances):
    #         if distance <= radius:
    #             al = self.atom_list_by_kind()
    #             neighbors.append(self.atom_list_by_kind('*')[index])
    #         else:
    #             break
    #
    #     return neighbors

    # # @property
    # def coords_nparray(self):
    #     atoms = self.atom_list_by_kind('*', excludeG=True)
    #     coords = np.ndarray(shape=(len(atoms), 3), dtype='float')
    #
    #     for idx, atom in enumerate(atoms):
    #         coords[idx] = atom.pos
    #
    #     return coords
    #
    #
    #
    # @property
    # def types_nparray(self):
    #     atoms = self.atom_list_by_kind('*', excludeG=True)
    #     types = np.empty(len(atoms), dtype='string')
    #
    #     for idx, atom in enumerate(atoms):
    #         types[idx] = atom.kind
    #
    #     return types
    #
    # @property
    # def bonds_nparray(self):
    #     atoms, atom_id_to_idx = self.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)
    #     bonds = self.bond_list_by_kind('*')
    #
    #     bonds_nparray = np.ndarray(shape=(len(bonds), 2), dtype='int')
    #     bond_types = np.empty(len(bonds), dtype='string')
    #
    #     for idx, bond in enumerate(bonds):
    #         bonds_nparray[idx, 0] = atom_id_to_idx[id(bond._atom1)]
    #         bonds_nparray[idx, 1] = atom_id_to_idx[id(bond._atom2)]
    #         bond_types[idx] = bond.kind
    #
    #     # return bonds_nparray, bond_types
    #     return bonds_nparray


    # def boundingbox(self, excludeG=True):
    #     """Compute the bounding box of the compound.
    #
    #     Returns:
    #         Box: Simulation box initialzied with min and max coordinates.
    #
    #     """
    #     minx = np.inf
    #     miny = np.inf
    #     minz = np.inf
    #     maxx = -np.inf
    #     maxy = -np.inf
    #     maxz = -np.inf
    #
    #     for a in self.atoms():
    #         if excludeG and a.kind == 'G':
    #             continue
    #         if a.pos[0] < minx:
    #             minx = a.pos[0]
    #         if a.pos[0] > maxx:
    #             maxx = a.pos[0]
    #         if a.pos[1] < miny:
    #             miny = a.pos[1]
    #         if a.pos[1] > maxy:
    #             maxy = a.pos[1]
    #         if a.pos[2] < minz:
    #             minz = a.pos[2]
    #         if a.pos[2] > maxz:
    #             maxz = a.pos[2]
    #
    #     min_coords = np.array([minx, miny, minz])
    #     max_coords = np.array([maxx, maxy, maxz])
    #
    #     return Box(mins=min_coords, maxes=max_coords)

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        if len(memo) == 0:
            memo[0] = self
        memo[id(self)] = newone

        newone.kind = deepcopy(self.kind, memo)
        newone.periodicity = deepcopy(self.periodicity, memo)

        # Copy the parent of everybody, except topmost compound being deepcopied.
        if memo[0] == self:
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)

        # Copy parts, except bonds with atoms outside the hierarchy.
        newone.parts = OrderedSet()
        for part in self.parts:
            if isinstance(part, Bond):
                if memo[0] in part.atom1.ancestors() and memo[0] in part.atom2.ancestors():
                    newone.parts.add(deepcopy(part,memo))
            else:
                newone.parts.add(deepcopy(part,memo))

        # Copy labels, except bonds with atoms outside the hierarchy
        newone.labels = OrderedDict()
        for k, v in self.labels.items():
            if isinstance(v, Bond):
                if memo[0] in v.atom1.ancestors() and memo[0] in v.atom2.ancestors():
                    newone.labels[k] = deepcopy(v, memo)
                    newone.labels[k].referrers.add(newone)
            else:
                newone.labels[k] = deepcopy(v, memo)
                if not isinstance(newone.labels[k], list):
                    newone.labels[k].referrers.add(newone)

        # Copy referrers that do not point out of the hierarchy.
        newone.referrers = set()
        for r in self.referrers:
            if memo[0] in r.ancestors():
                newone.referrers.add(deepcopy(r,memo))

        return newone
