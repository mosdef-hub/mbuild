from mbuild.atom import Atom
from mbuild.bond import Bond
from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.periodic_kdtree import PeriodicCKDTree

__author__ = 'sallai'
import numpy as np

class System(object):

    @classmethod
    def from_compound(cls, compound):
        atom_list, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

        n_atoms = len(atom_list)
        coords = np.ndarray(shape=(n_atoms, 3), dtype='float')
        types = np.empty(n_atoms, dtype='string')

        for idx, atom in enumerate(atom_list):
            coords[idx] = atom.pos
            types[idx] = atom.kind

        bond_list, bond_id_to_idx = compound.bond_list_by_kind('*', with_id_to_idx_mapping=True)
        n_bonds = len(bond_list)

        bonds = np.ndarray(shape=(n_bonds, 2), dtype='int')
        bond_types = np.empty(len(bonds), dtype='string')

        for idx, bond in enumerate(bond_list):
            bonds[idx, 0] = atom_id_to_idx[id(bond._atom1)]
            bonds[idx, 1] = atom_id_to_idx[id(bond._atom2)]
            bond_types[idx] = bond.kind

        sys = cls(coords=coords, types=types, bonds=bonds, bond_types=bond_types)
        sys.atom_list = atom_list
        sys.bond_list = bond_list
        sys.atom_id_to_idx = atom_id_to_idx
        sys.bond_id_to_idx = bond_id_to_idx
        return sys


    def __init__(self, coords=None, masses=None, charges=None, types=None, bonds=None, bond_types=None, angles=None, dihedrals=None, impropers=None, box=None):
        if coords is not None:
            self.coords = np.asarray(coords, 'float')
            self.n_atoms = self.coords.shape[0]
            assert(self.coords.shape[1] == 3)
        else:
            self.coords = None
            self.n_atoms = 0

        if masses is not None:
            self.masses = np.asarray(masses, 'float')
            assert(len(self.masses.shape) == 1)
        else:
            self.masses = None

        if charges is not None:
            self.charges= np.asarray(charges, 'float')
            assert(len(self.charges.shape) == 1)
        else:
            self.charges = None

        if types is not None:
            self.types= np.asarray(types, 'str')
            assert(len(self.types.shape) == 1)
        else:
            self.types = None

        if bonds is not None:
            self.bonds= np.asarray(bonds, 'int')
            self.n_bonds = self.bonds.shape[0]
            assert(self.bonds.shape[1] == 2)
        else:
            self.bonds = None
            self.n_bonds = 0

        if bond_types is not None:
            self.bond_types= np.asarray(bond_types, 'string')
            assert(len(self.bond_types.shape) == 1)
        else:
            self.bond_types = None

        if angles is not None:
            self.angles= np.asarray(angles, 'int')
            assert(self.angles.shape[1] == 3)
        else:
            self.angles = None

        if dihedrals is not None:
            self.dihedrals= np.asarray(dihedrals, 'int')
            assert(self.dihedrals.shape[1] == 4)
        else:
            self.dihedrals = None

        if impropers is not None:
            self.impropers= np.asarray(impropers, 'int')
            assert(self.impropers.shape[1] == 4)
        else:
            self.impropers = None

        if box is not None:
            self.box = box
            assert(isinstance(self.box, Box))
        else:
            self.box = Box()


    def update_compound(self, compound):
        atoms, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

        assert (len(atoms) == self.n_atoms)

        for idx, atom in enumerate(atoms):
            atom.pos = self.coords[idx]

    def to_compound(self, part=None):
        if part is None:
            part = Compound()

        atom_list = []
        for idx,kind in enumerate(self.types):
            coord = self.coords[idx]
            new_atom = Atom(str(kind), coord)
            part.add(new_atom, label="{0}[$]".format(kind))
            part.add(new_atom, label="atom[$]", containment=False)
            atom_list.append(new_atom)

        if self.bonds is not None:
            for idx,bond in enumerate(self.bonds):
                atom1 = atom_list[bond[0]]
                atom2 = atom_list[bond[1]]
                if(self.bond_types):
                    kind = self.bond_types[idx]
                else:
                    kind = None
                part.add(Bond(atom1, atom2, kind=kind), label="bond[$]")

        return part

    def boundingbox(self):
        return Box(mins=np.amin(self.coords, axis=0), maxes=np.amax(self.coords, axis=0))

    def _init_atom_kdtree(self):
            if len(self.coords) > 0:
                self._atom_kdtree = PeriodicCKDTree(self.coords)
            else:
                self._atom_kdtree = None

    def atoms_in_range(self, point, radius, maxItems=10):

        # create kdtree if it's not yet there
        if not hasattr(self,'_atom_kdtree'):
            self._init_atom_kdtree()

        if self._atom_kdtree is None:
            return []

        distances, indices = self._atom_kdtree.query(point, maxItems)

        neighbors = []
        for index, distance in zip(indices, distances):
            if distance <= radius:
                neighbors.append(self.atom_list[index])
            else:
                break

        return neighbors
