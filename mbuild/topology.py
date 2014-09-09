from mdtraj.core.element import Element

__author__ = 'sallai'
from mdtraj.core.topology import Topology as MDTTopology
from mdtraj.core import element as elem

class Topology(MDTTopology):

    def __init__(self):
        super(Topology, self).__init__()
        self._angles = []
        self._dihedrals = []


    def add_angle(self, atom1, atom2, atom3):
        print "Adding angle: {}-{}-{}".format(atom1, atom2, atom3)
        self._angles.append((atom1, atom2, atom3))

    @classmethod
    def from_compound(cls, compound, atom_list=None, bond_list=None):

        out = cls()
        atom_mapping = {}

        c = out.add_chain()
        r = out.add_residue("RES", c)

        if atom_list is None:
            atom_list, atom_id_to_idx = compound.atom_list_by_kind('*', excludeG=True, with_id_to_idx_mapping=True)

        for atom in atom_list:
            e = None
            try:
                e = elem.get_by_symbol(atom.kind)
            except:
                e = Element(  1000, atom.kind, atom.kind, 1.0)

            a = out.add_atom(atom.kind, e, r)
            atom_mapping[atom] = a

        if bond_list is None:
            bond_list, bond_id_to_idx = compound.bond_list_by_kind('*', with_id_to_idx_mapping=True)
        for idx, bond in enumerate(bond_list):
            a1 = bond.atom1
            a2 = bond.atom2
            out.add_bond(atom_mapping[a1], atom_mapping[a2])

        return out
