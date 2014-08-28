from mbuild.atom import Atom
from mbuild.bond import Bond
from mbuild.box import Box
from mbuild.compound import Compound

__author__ = 'sallai'
import numpy as np

class System(object):

    def __init__(self, coords=None, masses=None, charges=None, types=None, bonds=None, angles=None, dihedrals=None, impropers=None, box=None, compound=None):
        if coords is not None:
            self.coords = np.asarray(coords, 'float')
            assert(self.coords.shape[1] == 3)
        else:
            self.coords = None

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
            assert(self.bonds.shape[1] == 2)
        else:
            self.bonds = None

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

        # initialize from compound
        if compound is not None:

            atoms, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

            self.n_atoms = len(atoms)
            self.coords = np.ndarray(shape=(self.n_atoms, 3), dtype='float')
            self.types = np.empty(self.n_atoms, dtype='string')

            for idx, atom in enumerate(atoms):
                self.coords[idx] = atom.pos
                self.types[idx] = atom.kind

            bonds = compound.bond_list_by_kind('*')
            self.n_bonds = len(bonds)

            self.bonds = np.ndarray(shape=(len(bonds), 2), dtype='int')
            self.bond_types = np.empty(len(bonds), dtype='string')

            for idx, bond in enumerate(bonds):
                self.bonds[idx, 0] = self.atom_id_to_idx[id(bond._atom1)]
                self.bonds[idx, 1] = self.atom_id_to_idx[id(bond._atom2)]
                self.bond_types[idx] = bond.kind


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
            print new_atom
            part.add(new_atom, label="{0}[$]".format(kind))
            part.add(new_atom, label="atom[$]", containment=False)
            atom_list.append(new_atom)

        if self.bonds is not None:
            for idx,bond in enumerate(self.bonds):
                atom1 = atom_list[bond[0]]
                atom2 = atom_list[bond[1]]
                kind = self.bond_types[idx]
                part.add(Bond(atom1, atom2), kind=kind)

        return part
