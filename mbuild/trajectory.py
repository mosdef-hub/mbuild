import sys
import os

import numpy as np
import mdtraj as md

from mbuild.atom import Atom
from mbuild.bond import Bond
from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.formats.hoomdxml import save_hoomdxml
from mbuild.periodic_kdtree import PeriodicCKDTree
from mbuild.topology import Topology
from mbuild.formats.gromacs import save_gromacs
from mbuild.formats.lammps_data import save_lammps_data
from mbuild.formats.mol2 import save_mol2
from mbuild.formats.xyz import save_xyz


class Trajectory(object):
    """Derivative of MDTraj's Trajectory class with additional functionalities.

    Most notably, provides conversion to and from mBuild Compounds.
    """
    def __init__(self, *args, **kwargs):
        self._atom_kdtrees = {}

        if "trajectory" in kwargs:
            # We're wrapping an mdtraj Trajectory passed as a kwarg.
            trajectory = kwargs["trajectory"]
            self._w_trajectory = trajectory
        else:
            # We're wrapping an newly created md.Trajectory.
            self._w_trajectory = md.Trajectory(*args, **kwargs)

        # Wrap the topology.
        self._w_trajectory.topology = Topology(topology=self._w_trajectory.topology)

    @property
    def topology(self):
        return self._w_trajectory.topology

    @topology.setter
    def topology(self, t):
        if hasattr(t, '_w_topology'):
            # It's a wrapper -- use it as is.
            self._w_trajectory.topology = t
        else:
            assert isinstance(t, md.Topology)
            # It's an mdtraj topology -- wrap it.
            self._w_trajectory.topology = Topology(topology=t)

    @property
    def top(self):
        return self._topology

    @top.setter
    def top(self, t):
        self.topology = t

    @classmethod
    def from_compound(cls, compound, show_ports=False):
        exclude = not show_ports
        atom_list = compound.atom_list_by_kind('*', excludeG=exclude)

        t = Topology.from_compound(compound, atom_list=atom_list)

        n_atoms = len(atom_list)
        xyz = np.ndarray(shape=(1, n_atoms, 3), dtype='float')

        for idx, atom in enumerate(atom_list):
            xyz[0, idx] = atom.pos

        box = compound.boundingbox()
        unitcell_lengths = np.empty(3)
        for dim, val in enumerate(compound.periodicity):
            if val:
                unitcell_lengths[dim] = val
            else:
                unitcell_lengths[dim] = box.lengths[dim]
        traj = cls(xyz, t, unitcell_lengths=unitcell_lengths,
                   unitcell_angles=np.array([90, 90, 90]))
        return traj

    def update_compound(self, compound, frame=0):
        assert(isinstance(compound, Compound))

        for chain in self.topology.chains:
            for res in chain.residues:
                for atom in res.atoms:
                    compound.atoms[atom.index].pos = self.xyz[frame, atom.index]

    def to_compound(self, part=None, frame=0):
        if part is None:
            part = Compound()

        assert(isinstance(part, Compound))
        atom_mapping = {}
        for chain in self.topology.chains:
            if self.topology.n_chains > 1:
                chain_compound = Compound()
                part.add(chain_compound, "chain[$]")
            else:
                chain_compound = part
            for res in chain.residues:
                for atom in res.atoms:
                    new_atom = Atom(str(atom.name), self.xyz[frame, atom.index])
                    chain_compound.add(new_atom, label="{0}[$]".format(atom.name))
                    atom_mapping[atom] = new_atom

        for a1, a2 in self.topology.bonds:
            atom1 = atom_mapping[a1]
            atom2 = atom_mapping[a2]
            part.add(Bond(atom1, atom2))

        if (np.any(self.unitcell_lengths) and np.any(self.unitcell_lengths[0])):
            part.periodicity = self.unitcell_lengths[0]
        else: 
            part.periodicity = np.array([0., 0., 0.])
        return part

    def append_compound(self, compound):
        assert(isinstance(compound, Compound))
        self.to_compound(part=compound)

    def boundingbox(self, step):
        """Find the Box that encloses the all Atoms in a frame. """
        mins = np.amin(self.xyz[step], axis=0)
        maxes = np.amax(self.xyz[step], axis=0)
        return Box(mins=mins, maxs=maxes)

    def _init_atom_kdtree(self, frame=0):
            if self.n_atoms > 0:
                self._atom_kdtrees[frame] = PeriodicCKDTree(self.xyz[frame])
            else:
                self._atom_kdtrees[frame] = None

    def atoms_in_range_idx(self, point, radius, max_items=10, frame=0):
        """Return the indices of Atoms within a radius of a point.

        Args:
            point (list): The reference point in cartesian coordinates.
            radius (float): Find Atoms within this radius.
            max_items (int): Maximum number of Atoms to find.
            frame (int): The frame in the trajectory to operate on.
        Returns:
            neighbors (list): Indices of Atoms within specified range.
        """
        # Create kdtree if it's not yet there.
        if not frame in self._atom_kdtrees:
            self._init_atom_kdtree(frame=frame)
        if self._atom_kdtrees[frame] is None:
            return []

        distances, indices = self._atom_kdtrees[frame].query(point, max_items)
        neighbors = []
        for index, distance in zip(indices, distances):
            if distance <= radius:
                neighbors.append(index)
            else:
                # `distances` is sorted so we can exit the loop here.
                break
        return neighbors

    def atoms_in_range_xyz(self, *args, **kwargs):
        """Return coordinates of Atoms within a radius of a point. """
        neighbors_idx = self.atoms_in_range_idx(*args, **kwargs)
        return [self.xyz[kwargs["frame"], idx] for idx in neighbors_idx]

    def atoms_in_range(self, *args, **kwargs):
        """Return Atoms within a radius of a point. """
        neighbors_idx = self.atoms_in_range_idx(*args, **kwargs)
        return [self.topology.atom(idx) for idx in neighbors_idx]

    def bonds_by_atom_type(self, type_A, type_B):
        """Return a list of Bonds that contain Atoms of type_A and type_B. """
        bond_list = []
        for ab1 in self.top.bonds:
            ab = None
            if ab1[0].name == type_A and ab1[1].name == type_B:
                ab = ab1
            if ab1[1].name == type_A and ab1[0].name == type_B:
                ab = ab1[::-1]
            if ab is None:
                continue
            bond_list.append(ab)
        return bond_list

    @classmethod
    def load(cls, filename, relative_to_module=None, **kwargs):
        """Load a Trajectory from a file.

        Args:
            filename (str): Name of the input file.
            relative_to_module (bool, optional):  Look for data file in same
                directory as this python module.
        Returns:
            Trajectory
        """
        if relative_to_module is not None:
            current_dir = os.path.dirname(os.path.realpath(sys.modules[relative_to_module].__file__))
            filename = os.path.join(current_dir, filename)

        t = md.load(filename, **kwargs)
        return cls(trajectory=t, **kwargs)

    def save(self, filename, **kwargs):
        """Save the Trajectory to a file.

        Adds additional file formats beyond those provided by MDTraj.

        TODO: Write proper file classes that can be submitted to MDTraj.
        Args:
            filename (str): Name of the output file.
        """
        if filename.endswith(".hoomdxml"):
            save_hoomdxml(traj=self, filename=filename, **kwargs)
        elif filename.endswith(".gro") or filename.endswith(".top"):
            basename = ''.join(filename.split('.')[:-1])
            save_gromacs(traj=self, basename=basename, **kwargs)
        elif filename.endswith(".xyz"):
            save_xyz(traj=self, filename=filename, **kwargs)
        elif filename.endswith(".mol2"):
            save_mol2(traj=self, filename=filename, **kwargs)
        elif filename.startswith("data.") or filename.startswith(".lmp"):
            save_lammps_data(traj=self, filename=filename, **kwargs)
        else:
            self._w_trajectory.save(filename, **kwargs)

    def __getattr__(self, attr_name):
        return getattr(self._w_trajectory, attr_name)

    def __setattr__(self, key, value):
        if key in ['unitcell_vectors', 'unitcell_lengths', 'unitcell_angles', 'xyz', 'time']:
            self._w_trajectory.__setattr__(key, value)
        self.__dict__[key] = value

    def __getitem__(self, key):
        """Get a slice of this trajectory"""
        return self._w_trajectory.slice(key)

    def __str__(self):
        return "<%s>" % (self._string_summary_basic())

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))

    def _string_summary_basic(self):
        """Basic summary of traj in string form."""
        unitcell_str = 'and unitcells' if self._have_unitcell else 'without unitcells'
        value = "mbuild.Trajectory with %d frames, %d atoms, %d residues, %s" % (
                    self.n_frames, self.n_atoms, self.n_residues, unitcell_str)
        return value
