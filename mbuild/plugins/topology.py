__author__ = 'sallai'
from mdtraj.core.topology import Topology as MDTTopology
from mdtraj.core import element as elem

class Topology(MDTTopology):

    @classmethod
    def from_compound(cls, compound, atom_list=None, bond_list=None):

        out = cls()
        atom_mapping = {}

        c = out.add_chain()
        r = out.add_residue("RES", c)

        if atom_list is None:
            atom_list, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

        for atom in atom_list:
            a = out.add_atom(atom.kind, elem.get_by_symbol(atom.kind), r)
            atom_mapping[atom] = a

        if bond_list is None:
            bond_list, bond_id_to_idx = compound.bond_list_by_kind('*', with_id_to_idx_mapping=True)
        for idx, bond in enumerate(bond_list):
            a1 = bond.atom1
            a2 = bond.atom2
            out.add_bond(atom_mapping[a1], atom_mapping[a2])

        return out
