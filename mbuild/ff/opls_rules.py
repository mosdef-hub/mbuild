__author__ = 'sallai'

from itertools import combinations, product
import pdb

from mbuild.plot import Plot
from mbuild.prototype import Prototype
from mbuild.rules import RuleEngine
import mbuild.unit as units
from wip.mpc_monolayer.mpc_monomer import MpcMonomer
from opls_forcefield import OplsForceField


class OplsRules(RuleEngine):
    """Loads a forcefield and parameterizes a compound."""
    def __init__(self, compound, force_field):
        super(OplsRules, self).__init__(compound)
        self.force_field = force_field

    def execute(self):
        r_err = 0.4 * units.angstroms
        unique_bond_types = set()

        for atom in self.compound.atoms():
            if atom.kind != "G":
                unique_bond_types.add(Prototype.getAttr(atom.kind,
                        "bond_type"))
                try:
                    unique_bond_types.add(Prototype.getAttr(atom.kind,
                        "bond_type")[0]+"*")
                except:
                    pdb.set_trace()


        pairs = [(i, j) for i,j in combinations(unique_bond_types, 2)]

        for pair in pairs:
            if pair in self.force_field.bond_types:
                r, k = self.force_field.bond_types[pair]
                print "Found for pair {0}".format(pair)
            elif pair[::-1] in self.force_field.bond_types:
                r, k = self.force_field.bond_types[pair[::-1]]
                print "Found for pair {0}".format(pair)
            else:
                print "No bond type found for pair {0}".format(pair)
                continue

            # pair[0] and pair[1] are bond type labels

            # find all atoms that match the bond type labels (kind, alias with optional wildcards)
            atomTypes1 = self.force_field.findAtomTypes(pair[0])
            atomTypes2 = self.force_field.findAtomTypes(pair[1])


            # for every combination of the matching atom kinds, create a bond type
            for (atomType1, atomType2) in product(atomTypes1, atomTypes2):

                # print "Adding rule " + atomType1 + "-" + atomType2

                if(atomType1 < atomType2):
                    pair = (atomType1, atomType2)
                else:
                    pair = (atomType1, atomType2)

                # Get all atoms that form this pair and try to add a bond
                self.add_bond(pair[0], pair[1],
                        (r - r_err)._value, (r + r_err)._value,
                        "{0}-{1}".format(pair[0], pair[1]))

if __name__ == "__main__":
    m = MpcMonomer()
    ff = OplsForceField()
    ff.get_atom_types(m)

    rules = OplsRules(m, ff)
    rules.execute()
    #pdb.set_trace()
    Plot(m).show()


