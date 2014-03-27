from itertools import combinations, product
from scimath.units.mass import gram
from scimath.units.quantity import Quantity
from traits.trait_types import Instance, Dict, Tuple
from atomtype import AtomType
from bondtype import BondType
from mbuild.plot import Plot
from mbuild.prototype import Prototype
from forcefieldrules import ForceFieldRules
from wip.mpc_monolayer.mpc_monomer import MpcMonomer

__author__ = 'sallai'
import os
import pdb
from traits.api import HasTraits, Str, Float
from scimath.units.quantity_traits import QuantityTrait



class GromacsFF(HasTraits):
    """A reader for GROMACS forcefield files
    """

    atom_types = Dict(key_trait=Str, value_trait=Instance(AtomType))
    unresolved_bond_types = Dict(key_trait=Tuple, value_trait=Instance(BondType))
    bond_types = Dict(key_trait=Tuple, value_trait=Instance(BondType))

    def __init__(self, forcefield_dir):

        gmx_lib = os.getenv('GMXLIB', 'gromacs_forcefields')

        # read the nonbonded file
        atom_type_file = os.path.join(gmx_lib, forcefield_dir, 'ffnonbonded.itp')
        with open(atom_type_file, 'r') as f:
            for line in f:
                if line.strip():
                    line = line.strip()

                    if line[0] in [';', '#']:
                        # it's a comment or a preprocessor directive
                        continue

                    if line[0] in ['[']:
                        # its a section section_keyword
                        section_keyword = line
                        continue

                    if section_keyword == '[ atomtypes ]':
                        fields = line.split()
                        atomType = AtomType(kind=fields[0], alias=fields[1], atomicNumber=fields[2], mass=fields[3], charge=fields[4], sigma=fields[6], epsilon=fields[7])
                        self.atom_types[fields[0]] = atomType
                    else:
                        # print "Not parsing line in section " + section_keyword +":"+line
                        pass


        # read the bonded file
        bond_type_file = os.path.join(gmx_lib, forcefield_dir, 'ffbonded.itp')
        with open(bond_type_file, 'r') as f:
            for line in f:
                if line.strip():
                    line = line.strip()

                    if line[0] in [';', '#']:
                        # it's a comment or a preprocessor directive
                        continue

                    if line[0] in ['[']:
                        # its a section section_keyword
                        section_keyword = line
                        continue

                    if section_keyword == '[ bondtypes ]':
                        fields = line.split()
                        atomKindIdentifier1 = fields[0]
                        atomKindIdentifier2 = fields[1]

                        if(atomKindIdentifier1 < atomKindIdentifier2):
                            pair = (atomKindIdentifier1, atomKindIdentifier2)
                        else:
                            pair = (atomKindIdentifier2, atomKindIdentifier1)

                        self.unresolved_bond_types[pair] = BondType(bond_types=pair, r=float(fields[3])/10, k=fields[4])
                    else:
                        # print "Not parsing line in section " + section_keyword +":"+line
                        pass

        # resolve aliases and wildcards in bond types
        self.resolveBondTypes()

    def resolveBondTypes(self):
        for (atomKindIdentifier1, atomKindIdentifier2), bondType in self.unresolved_bond_types.iteritems():

            # find all atoms that match the atomKindIdentifiers (kind, alias with optional wildcards)
            atomTypes1 = self.findAtomTypes(atomKindIdentifier1)
            atomTypes2 = self.findAtomTypes(atomKindIdentifier2)

            # for every combination of the matching atom kinds, create a bond type
            for (atomType1, atomType2) in product(atomTypes1, atomTypes2):

                if(atomType1.kind < atomType2.kind):
                    pair = (atomType1.kind, atomType2.kind)
                else:
                    pair = (atomType2.kind, atomType1.kind)

                self.bond_types[pair] = BondType(bond_types=pair, r=bondType.r, k=bondType.k)


    def findAtomTypes(self, id):
        # if the id is the atom type, return the AtomType object
        if id in self.atom_types:
            return [self.atom_types[id]]

        matching_atom_types = []

        for at in self.atom_types.values():
            if str(id).endswith('*'):
                # id is a wildcard ending in *
                prefix = str(id)[0:-1]

                if at.alias.startswith(prefix):
                    matching_atom_types.append(at)
                elif at.kind.startswith(prefix):
                    matching_atom_types.append(at)
            else:
                # id is not a wildcard
                if id == at.alias:
                    matching_atom_types.append(at)
                elif id == at.kind:
                    matching_atom_types.append(at)

        # return the matching_atom_types
        return matching_atom_types

    def rules(self):
        """
        Create a rule engine using the force field information
        """
        return ForceFieldRules(self)

if __name__ == "__main__":
    # for bt in ff.bond_types:
    #     print bt
    #     print ff.bond_types[bt].print_traits()

    m = MpcMonomer()
    ff = GromacsFF("oplsaa.ff")

    ff.rules().apply(m)

    print m.bonds

    Plot(m).show()


