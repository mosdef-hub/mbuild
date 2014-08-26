from mbuild.box import Box
from mbuild.compound import Compound

__author__ = 'sallai'
import numpy as np

class System(object):

    def __init__(self, coords=None, masses=None, charges=None, types=None, bonds=None, angles=None, dihedrals=None, impropers=None, box=None):
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


    def from_compound(self, compound):
        assert(isinstance(compound, Compound))

        n_atoms = len([atom for atom in compound.atoms() if atom.kind != "G"])
        n_bonds = len(list(compound.bonds()))

        self.coords = np.array(shape=(n_atoms,3))
        self.types = np.array(shape=(n_atoms,1))

        id_to_idx = dict()
        for atom_idx, atom in enumerate([a for a in compound.atoms() if a.kind != "G"]):
            id_to_idx[id(atom)] = atom_idx
            self.coords[atom_idx] = atom.pos
            self.types[atom_idx] = atom.kind

        if n_bonds:
            for bond_idx, bond in enumerate(compound.bonds()):
                self.bonds[bond_idx][0] = id_to_idx[id(bond.atom1)]
                self.bonds[bond_idx][1] = id_to_idx[id(bond.atom2)]
