import pdb
from itertools import product
from mbuild.prototype import Prototype

__author__ = 'sallai'

class ForceField(object):

    def __init__(self, name):
        """
        """
        self.name = name
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
        ff = ForceField(self.name)
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
        for (alias1, alias2), bond_type in self.bond_types.iteritems():
            # find all atoms that match the atomKindIdentifiers (kind, alias with optional wildcards)
            atom_types1 = ff.find_atom_types(alias1)
            atom_types2 = ff.find_atom_types(alias2)

            # for every combination of the matching atom kinds, create a bond type
            for (atom_type1, atom_type2) in product(atom_types1, atom_types2):

                # Why were we sorting these? Pretty sure it's unneccesary at this point
                # Also, this won't hold up for other forcefields
                #if atom_type1 < atom_type2:
                    #pair = (atom_type1, atom_type2)
                #else:
                #    pair = (atom_type2, atom_type1)

                pair = (atom_type1, atom_type2)
                ff.bond_types[pair] = bond_type

        # prune the angle types, resolving wildcards
        for (alias1, alias2, alias3), angle_type in self.angle_types.iteritems():
            # find all atoms that match the atomKindIdentifiers (kind, alias with optional wildcards)
            atom_types1 = ff.find_atom_types(alias1)
            atom_types2 = ff.find_atom_types(alias2)
            atom_types3 = ff.find_atom_types(alias3)

            # for every combination of the matching atom kinds, create a bond type
            for (atom_type1, atom_type2, atom_type3) in product(atom_types1, atom_types2, atom_types3):
                triplet = (atom_type1, atom_type2, atom_type3)
                #if atom_type1 < atom_type2:
                #    triplet = (atom_type1, atom_type2)
                #else:
                #    triplet = (atom_type2, atom_type1)
                ff.angle_types[triplet] = angle_type

        # prune the dihedral types, resolving wildcards
        for (alias1, alias2, alias3, alias4), dihedral_type in self.dihedral_types.iteritems():
            # find all atoms that match the atomKindIdentifiers (kind, alias with optional wildcards)
            atom_types1 = ff.find_atom_types(alias1)
            atom_types2 = ff.find_atom_types(alias2)
            atom_types3 = ff.find_atom_types(alias3)
            atom_types4 = ff.find_atom_types(alias4)

            # for every combination of the matching atom kinds, create a bond type
            for (atom_type1, atom_type2, atom_type3, atom_type4) in product(atom_types1, atom_types2, atom_types3, atom_types4):
                quartet = (atom_type1, atom_type2, atom_type3, atom_type4)
                #if atom_type1 < atom_type2:
                #    triplet = (atom_type1, atom_type2)
                #else:
                #    triplet = (atom_type2, atom_type1)
                ff.dihedral_types[quartet] = dihedral_type

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


