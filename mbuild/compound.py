from collections import OrderedDict
from copy import deepcopy

import numpy as np

from atom import Atom
from bond import Bond
from box import Box
from mbuild.mbase import MBase
from mbuild.has_parts_mixin import HasPartsMixin
from mbuild.part_mixin import PartMixin
from periodic_kdtree import PeriodicCKDTree
from orderedset import OrderedSet

class Compound(MBase, PartMixin, HasPartsMixin):
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

        super(Compound, self).__init__()

        # Set kind to classname if not specified.
        if kind:
            self.kind = kind
        else:
            self.kind = self.__class__.__name__

        # A periodocity of zero in any direction is treated as non-periodic.
        if not periodicity:
            periodicity = np.array([0.0, 0.0, 0.0])
        self.periodicity = periodicity

    def atoms(self):
        return self._yield_parts(Atom)

    def bonds(self):
        return self._yield_parts(Bond)


    def referenced_ports(self):
        from mbuild.port import Port
        return [port for port in self.labels.values() if isinstance(port, Port)]


    def post_remove(self, removed_part):
        super(Compound, self).post_remove(removed_part)

        # If removing an atom, make sure to remove the bonds it's part of.
        if isinstance(removed_part, Atom):
            for bond in removed_part.bonds:
                if bond.parent is not None:
                    bond.parent.remove(bond)


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


    def append_from_file(self, filename, relative_to_module=None, frame=0):

        # Load trajectory from file
        from plugins.trajectory import Trajectory
        traj = Trajectory.load(filename, relative_to_module=relative_to_module)

        self.append_from_trajectory(traj, frame=frame)

    def update_from_file(self, filename, relative_to_module=None, frame=0):

        # Load trajectory from file
        from plugins.trajectory import Trajectory
        traj = Trajectory.load(filename, relative_to_module=relative_to_module)

        self.update_from_trajectory(traj, frame=frame)

    def append_from_trajectory(self, traj, frame=0):
        # Append to this compound the trajectory's topology (atoms, bonds) with the atom's xyz positions as of frame 0
        traj.to_compound(part=self, frame=frame)

    def update_from_trajectory(self, traj, frame=0):
        traj.update_compound(self, frame=frame)

    def save(self, filename, **kwargs):
        self.to_trajectory().save(filename, **kwargs)

    def to_trajectory(self):
        from plugins.trajectory import Trajectory
        return Trajectory.from_compound(self)

    @classmethod
    def load(cls, filename, relative_to_module=None, frame=0):
        from plugins.trajectory import Trajectory
        traj = Trajectory.load(filename, relative_to_module=relative_to_module)
        return traj.to_compound(frame=frame)

    def min_periodic_distance(self, x0, x1):
        """Vectorized distance calculation considering minimum image. """
        d = np.abs(x0 - x1)
        d = np.where(d > 0.5 * self.periodicity, self.periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

    def boundingbox(self, excludeG=True):
        """Compute the bounding box of the compound.

        Returns:
            Box: Simulation box initialzied with min and max coordinates.

        """
        minx = np.inf
        miny = np.inf
        minz = np.inf
        maxx = -np.inf
        maxy = -np.inf
        maxz = -np.inf

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
