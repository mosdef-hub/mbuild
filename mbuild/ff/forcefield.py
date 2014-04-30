import pdb
from itertools import product
from mbuild.prototype import Prototype

__author__ = 'sallai'

class ForceField(object):

    def __init__(self):
        """

        """
        self.atom_types = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.dihedral_types = dict()

    def find_atom_types(self, atom_id):
        # if the id is the atom type, return the AtomType object
        """
        if atom_id == 'SI': pdb.set_trace()
        if atom_id in self.atom_types:
            return [self.atom_types[atom_id].kind]
        """
        matching_atom_types = []

        for kind, atom_type in self.atom_types.iteritems():

            if str(atom_id).endswith('*'):
                # id is a wildcard ending in *
                prefix = str(atom_id)[:-1]

                if atom_type.alias.startswith(prefix):
                    matching_atom_types.append(kind)
                elif kind.startswith(prefix):
                    matching_atom_types.append(kind)
            else:
                # id is not a wildcard
                if atom_id == atom_type.alias:
                    matching_atom_types.append(kind)

        # return the matching_atom_types
        return matching_atom_types

    def prune(self, compound, verbose=False):
        """
        Returns a pruned force field that contains information only relevant to the compound.
        """

        # create new force field object
        ff = ForceField()
        ff.atom_types.update(self.atom_types)

        all_kinds = ff.atom_types.keys()
        retained_kinds = compound.atomKinds()

        # prune the atom types
        for atomkind in all_kinds:
            if atomkind not in retained_kinds:
                del ff.atom_types[atomkind]
                if verbose:
                    print "removing " + atomkind
            else:
                if verbose:
                    print "retaining " + atomkind

        # prune the bond types, resolving wildcards
        for (alias1, alias2), bondType in self.bond_types.iteritems():
            # find all atoms that match the atomKindIdentifiers (kind, alias with optional wildcards)
            atomTypes1 = ff.find_atom_types(alias1)
            atomTypes2 = ff.find_atom_types(alias2)

            # for every combination of the matching atom kinds, create a bond type
            for (atomType1, atomType2) in product(atomTypes1, atomTypes2):

                if atomType1 < atomType2:
                    pair = (atomType1, atomType2)
                else:
                    pair = (atomType2, atomType1)
                ff.bond_types[pair] = bondType

        # set up prototypes
        ff.init_prototypes()

        return ff

    def init_prototypes(self):
        """Parameterize a compound with atom specific information."""
        for kind in self.atom_types:
            params = self.atom_types[kind]
            Prototype(kind,
                    bond_type=params.alias,
                    atomic_num=params.atomic_number,
                    mass=params.mass,
                    charge=params.charge,
                    sigma=params.sigma,
                    epsilon=params.epsilon)


