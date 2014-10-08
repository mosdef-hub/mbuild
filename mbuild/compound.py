from collections import OrderedDict
from copy import deepcopy

import numpy as np

from atom import Atom
from bond import Bond
from box import Box
from mbuild.coordinate_transform import translate
from mbuild.mbase import MBase
from mbuild.has_parts_mixin import HasPartsMixin
from mbuild.part_mixin import PartMixin
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

    def n_atoms(self):
        return sum([1 for _ in self.atoms()])

    def n_bonds(self):
        return sum([1 for _ in self.bonds()])

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
        """Append to Compound with information from a Trajectory file. """
        from mbuild.trajectory import Trajectory
        traj = Trajectory.load(filename, relative_to_module=relative_to_module)

        self.append_from_trajectory(traj, frame=frame)

    def update_from_file(self, filename, relative_to_module=None, frame=0):
        """Update Compound with information from a Trajectory file. """
        from mbuild.trajectory import Trajectory
        traj = Trajectory.load(filename, relative_to_module=relative_to_module)

        self.update_from_trajectory(traj, frame=frame)

    def append_from_trajectory(self, traj, frame=0):
        """Append the Trajectory's topology (atoms, bonds). """
        traj.to_compound(part=self, frame=frame)

    def update_from_trajectory(self, traj, frame=0):
        """Update the Compound with the Trajectory's topology. """
        traj.update_compound(self, frame=frame)

    def save(self, filename, **kwargs):
        """Save the Compound to a file.

        This creates an intermediate Trajectory object.
        """
        self.to_trajectory().save(filename, **kwargs)

    def to_trajectory(self):
        """Convert the Compound to a Trajectory. """
        from mbuild.trajectory import Trajectory
        return Trajectory.from_compound(self)

    @classmethod
    def load(cls, filename, relative_to_module=None, frame=0):
        from mbuild.trajectory import Trajectory
        traj = Trajectory.load(filename, relative_to_module=relative_to_module)
        return traj.to_compound(frame=frame)

    def to_molecule(self):
        from pybel import Molecule
        from openbabel import OBMol
        from mbuild.prototype import Prototype

        mol = OBMol()

        atoms, atom_id_to_index = self.atom_list_by_kind(excludeG=True, with_id_to_idx_mapping=True)

        for atom in atoms:
            a = mol.NewAtom()

            atomic_num = Prototype.getAttr(atom.kind, "atomic_number", 0)
            print atomic_num

            a.SetAtomicNum(atomic_num)
            a.SetVector(float(atom.pos[0])*10, float(atom.pos[1])*10, float(atom.pos[2])*10) # coordinates

        for bond in self.bond_list_by_kind():
            idx1 = atom_id_to_index[id(bond.atom1)]
            idx2 = atom_id_to_index[id(bond.atom2)]
            mol.AddBond(idx1+1, idx2+1, 1)   # atoms indexed from 1

        return Molecule(mol)

    def update_from_molecule(self, mol):
        """

        Args:
            mol:
        """
        from pybel import Molecule
        assert(isinstance(mol, Molecule))

        atoms, atom_id_to_idx = self.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

        assert (len(atoms) == mol.OBMol.NumAtoms())

        idx = 0
        for atom in mol.atoms:
            print atom
            atoms[idx].pos = np.array(atom.coords) / 10.0
            idx += 1

    def visualize(self):
        """Visualize the Compound using VMD.

        Assumes you have VMD installed and can call it from the command line via
        'vmd'.

        TODO: Make more useful/robust. Look into pizza.py's vmd.py.
        """
        filename = 'visualize_{}.pdb'.format(self.__class__.__name__)
        traj = self.to_trajectory()
        traj.save(filename=filename)
        import os
        try:
            os.system('vmd {}'.format(filename))
        except OSError:
            print("Visualization with VMD failed. Make sure you it is installed"
                  "correctly and launchable from the command line via 'vmd'.")


    def min_periodic_distance(self, x0, x1):
        """Vectorized distance calculation considering minimum image. """
        d = np.abs(x0 - x1)
        d = np.where(d > 0.5 * self.periodicity, self.periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

    def boundingbox(self, excludeG=True):
        """Compute the bounding box of the compound.

        Args:
            excludeG (bool): Exclude Atoms of kind 'G' (typically reserved for
                             Ports).
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

        return Box(mins=min_coords, maxs=max_coords)

    def atoms_in_range(self, point, radius, max_items=10):
        """Find all Atoms within a radius of a point.

        Args:
            point:
            radius:
            max_items:
        Returns:
            list of Atoms within range
        """
        atoms = self.atom_list_by_kind(excludeG=True)
        traj = self.to_trajectory()
        idxs = traj.atoms_in_range_idx(point, radius, max_items=max_items)
        return [atoms[idx] for idx in idxs]

    def wrap(self):
        """ """
        assert np.any(self.periodicity)
        box = self.boundingbox()
        translate(self, -box.mins)
        for atom in self.atoms():
            for k, c in enumerate(atom.pos):
                if self.periodicity[k]:
                    if c < 0.0:
                        atom.pos[k] = self.periodicity[k] + c
                    if c > self.periodicity[k]:
                        atom.pos[k] = c - self.periodicity[k]


    def add_bond(self, type_a, type_b, dmin, dmax, kind=None):
        """ai-bj distance is in [dmin, dmax] => add bond a1xb(ai,bj) (symmetric)."""
        for a1 in self.atom_list_by_kind(type_a):
            nearest = self.atoms_in_range(a1.pos, dmax)
            for a2 in nearest:
                if (a2.kind==type_b) and (dmin <= self.min_periodic_distance(a2.pos, a1.pos) <= dmax):
                    self.add(Bond(a1, a2, kind=kind))

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        if len(memo) == 0:
            memo[0] = self
        memo[id(self)] = newone

        # first copy those attributes that don't need deepcopying
        newone.kind = deepcopy(self.kind, memo)
        newone.periodicity = deepcopy(self.periodicity, memo)

        # create empty containers
        newone.parts = OrderedSet()
        newone.labels = OrderedDict()
        newone.referrers = set()

        # Copy the parent of everybody, except topmost compound being deepcopied.
        if memo[0] == self:
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)

        # Copy parts, except bonds with atoms outside the hierarchy.
        for part in self.parts:
            if isinstance(part, Bond):
                if memo[0] in part.atom1.ancestors() and memo[0] in part.atom2.ancestors():
                    newone.parts.add(deepcopy(part,memo))
            else:
                newone.parts.add(deepcopy(part,memo))

        # Copy labels, except bonds with atoms outside the hierarchy
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
        for r in self.referrers:
            if memo[0] in r.ancestors():
                newone.referrers.add(deepcopy(r,memo))

        return newone
