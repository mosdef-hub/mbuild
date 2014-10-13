import sys
import os

import numpy as np
import mdtraj as md
from mdtraj.utils.six import string_types
from mdtraj.core.trajectory import load_frame, load_mol2, load_prmtop

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


class Trajectory(md.Trajectory):

    def __init__(self, *args, **kwargs):
        self._atom_kdtrees = {}

        if "trajectory" in kwargs:
            # we're casting an md.Trajectory to mbuild Trajectory
            trajectory = kwargs["trajectory"]
            assert(isinstance(trajectory, md.Trajectory))
            assert(len(kwargs) == 1)
            super(Trajectory, self).__init__(trajectory.xyz, trajectory.topology, time=trajectory.time, unitcell_lengths=trajectory.unitcell_lengths, unitcell_angles=trajectory.unitcell_angles)
        else:
            super(Trajectory, self).__init__(*args, **kwargs)

    @classmethod
    def from_compound(cls, compound):
        atom_list, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

        t = Topology.from_compound(compound, atom_list=atom_list)

        n_atoms = len(atom_list)
        xyz = np.ndarray(shape=(1, n_atoms, 3), dtype='float')

        for idx, atom in enumerate(atom_list):
            xyz[0,idx] = atom.pos

        box = compound.boundingbox()
        unitcell_lengths = np.empty(3)
        for dim, val in enumerate(compound.periodicity):
            if val:
                unitcell_lengths[dim] = val
            else:
                unitcell_lengths[dim] = box.lengths[dim]
        traj = cls(xyz, t, unitcell_lengths=unitcell_lengths, unitcell_angles=np.array([90, 90, 90]))
        return traj

    def update_compound(self, compound, frame=0):
        assert(isinstance(compound, Compound))

        atoms, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

        idx = 0
        for chain in self.topology.chains:
            for res in chain.residues:
                for atom in res.atoms:
                    atoms[idx].pos = self.xyz[frame, idx]
                    idx += 1


    def to_compound(self, part=None, frame=0):
        if part is None:
            part = Compound()

        assert(isinstance(part, Compound))

        atom_mapping = {}
        idx = 0
        for chain in self.topology.chains:
            if self.topology.n_chains > 1:
                chain_compound = Compound()
                part.add(chain_compound, "chain[$]")
            else:
                chain_compound = part
            for res in chain.residues:
                for atom in res.atoms:
                    new_atom = Atom(str(atom.name), self.xyz[frame, idx])
                    chain_compound.add(new_atom, label="{0}[$]".format(atom.name))
                    chain_compound.add(new_atom, label="atom[$]", containment=False)
                    atom_mapping[atom] = new_atom
                    idx += 1

        for a1, a2 in self.topology.bonds:
            atom1 = atom_mapping[a1]
            atom2 = atom_mapping[a2]
            part.add(Bond(atom1, atom2), label="bond[$]")

        return part

    def append_compound(self, compound):
        assert(isinstance(compound, Compound))

        self.to_compound(part=compound)

    def boundingbox(self, step):
        mins = np.amin(self.xyz[step], axis=0)
        maxes = np.amax(self.xyz[step], axis=0)
        return Box(mins=mins, maxs=maxes)

    def _init_atom_kdtree(self, frame=0):
            if self.n_atoms > 0:
                self._atom_kdtrees[frame] = PeriodicCKDTree(self.xyz[frame])
            else:
                self._atom_kdtrees[frame] = None

    def atoms_in_range_idx(self, point, radius, max_items=10, frame=0):

        # create kdtree if it's not yet there
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
                break

        return neighbors

    def atoms_in_range_xyz(self, *args, **kwargs):
        neighbors_idx = self.atoms_in_range_idx(*args, **kwargs)
        return [self.xyz[kwargs["frame"], idx] for idx in neighbors_idx]

    def atoms_in_range(self, *args, **kwargs):
        neighbors_idx = self.atoms_in_range_idx(*args, **kwargs)
        return [self.topology.atom(idx) for idx in neighbors_idx]

    def bonds_by_atom_type(self, type_A, type_B):
        bond_list = []
        for ab1 in self.top.bonds:
            ab = None
            if (ab1[0].name == type_A and ab1[1].name == type_B):
                ab  = ab1

            if (ab1[1].name == type_A and ab1[0].name == type_B):
                ab = ab1[::-1]

            if ab is None:
                continue

            bond_list.append(ab)

        return bond_list

    @staticmethod
    def order_bond(ab, type_A, type_B):
        """
        Order the atoms in a bond by atom type.
        :param ab: the bond (tuple)
        :param type_A: the type of the first atom
        :param type_B: the type of the second atom
        :return: ab with its atoms reordered such that its first atom is of type_A and its second atom is of type_B,
        or None if ab the atom types do not match
        """
        if ab[0].name == type_A and ab[1].name == type_B:
            return ab

        if ab[1].name == type_A and ab[0].name == type_B:
            return ab[::-1]
        return None

    def bonds_by_atom(self, atom):
        bond_list = []
        for bond in self.top.bonds:
            if atom in bond:
                bond_list.append(bond)
        return bond_list

    def neighbor_bonds(self, bond):
        atom0 = bond[0]
        atom1 = bond[1]
        atom0_bonds = set(self.bonds_by_atom(atom0))
        atom1_bonds = set(self.bonds_by_atom(atom1))
        return (atom0_bonds | atom1_bonds) - set([bond])

    @classmethod
    def load(cls, filename, relative_to_module=None, **kwargs):
        """

        Args:
            filename (str):
            relative_to_module (bool, optional):  Look for data file in same
                directory as this python module.
        Returns:
            Trajectory
        """
        if relative_to_module is not None:
            current_dir = os.path.dirname(os.path.realpath(sys.modules[relative_to_module].__file__))
            filename = os.path.join(current_dir, filename)

        t = md.load(filename)[0]  # TODO: figure out why this is a tuple
        return cls(trajectory=t, **kwargs)

    def save(self, filename, **kwargs):
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
            super(Trajectory, self).save(filename, **kwargs)


if __name__ == "__main__":
    t1 = Trajectory.load("../../../mbuild/tests/methyl.pdb")

    compound = t1.to_compound()

    print compound

    t2 = Trajectory.from_compound(compound)
    t2.save("../../../mbuild/tests/methyl2.pdb")
