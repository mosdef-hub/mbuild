from mdtraj.core.topology import Topology as MDTTopology
from mdtraj.core.trajectory import Trajectory as MDTTrajectory
from mbuild.atom import Atom
from mbuild.bond import Bond
from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.periodic_kdtree import PeriodicCKDTree
from mbuild.plugins.topology import Topology
import mdtraj as md

__author__ = 'sallai'
import numpy as np

class Trajectory(MDTTrajectory):

    def __init__(self, *args, **kwargs):
        if "trajectory" in kwargs:
            trajectory = kwargs["trajectory"]
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

        return cls(xyz, t)


    # def __init__(self, coords=None, masses=None, charges=None, types=None, bonds=None, bond_types=None, angles=None, dihedrals=None, impropers=None, box=None):
    #     if coords is not None:
    #         self.coords = np.asarray(coords, 'float')
    #         self.n_atoms = self.coords.shape[0]
    #         assert(self.coords.shape[1] == 3)
    #     else:
    #         self.coords = None
    #         self.n_atoms = 0
    #
    #     if masses is not None:
    #         self.masses = np.asarray(masses, 'float')
    #         assert(len(self.masses.shape) == 1)
    #     else:
    #         self.masses = None
    #
    #     if charges is not None:
    #         self.charges= np.asarray(charges, 'float')
    #         assert(len(self.charges.shape) == 1)
    #     else:
    #         self.charges = None
    #
    #     if types is not None:
    #         self.types= np.asarray(types, 'str')
    #         assert(len(self.types.shape) == 1)
    #     else:
    #         self.types = None
    #
    #     if bonds is not None:
    #         self.bonds= np.asarray(bonds, 'int')
    #         self.n_bonds = self.bonds.shape[0]
    #         assert(self.bonds.shape[1] == 2)
    #     else:
    #         self.bonds = None
    #         self.n_bonds = 0
    #
    #     if bond_types is not None:
    #         self.bond_types= np.asarray(bond_types, 'string')
    #         assert(len(self.bond_types.shape) == 1)
    #     else:
    #         self.bond_types = None
    #
    #     if angles is not None:
    #         self.angles= np.asarray(angles, 'int')
    #         assert(self.angles.shape[1] == 3)
    #     else:
    #         self.angles = None
    #
    #     if dihedrals is not None:
    #         self.dihedrals= np.asarray(dihedrals, 'int')
    #         assert(self.dihedrals.shape[1] == 4)
    #     else:
    #         self.dihedrals = None
    #
    #     if impropers is not None:
    #         self.impropers= np.asarray(impropers, 'int')
    #         assert(self.impropers.shape[1] == 4)
    #     else:
    #         self.impropers = None
    #
    #     if box is not None:
    #         self.box = box
    #         assert(isinstance(self.box, Box))
    #     else:
    #         self.box = Box()


    # def update_compound(self, compound):
    #     atoms, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)
    #
    #     assert (len(atoms) == self.n_atoms)
    #
    #     for idx, atom in enumerate(atoms):
    #         atom.pos = self.coords[idx]

    def to_compound(self, part=None, frame=0):
        if part is None:
            part = Compound()

        atom_mapping = {}
        idx = 0
        for chain in self.topology.chains:
            chain_compound = Compound()
            part.add(chain_compound, "chain[$]")
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

    # def boundingbox(self):
    #     return Box(mins=np.amin(self.coords, axis=0), maxes=np.amax(self.coords, axis=0))

    # def _init_atom_kdtree(self):
    #         if len(self.coords) > 0:
    #             self._atom_kdtree = PeriodicCKDTree(self.coords)
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
    #             neighbors.append(self.atom_list[index])
    #         else:
    #             break
    #
    #     return neighbors



def load(filename):
    t = md.load(filename)
    return Trajectory(trajectory=t)


if __name__ == "__main__":
    import mdtraj as md
    t1 = load("../../../mbuild/tests/methyl.pdb")

    compound = t1.to_compound()

    print compound

    t2 = Trajectory.from_compound(compound)
    t2.save("../../../mbuild/tests/methyl2.pdb")
