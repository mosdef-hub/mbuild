from itertools import combinations, product
from warnings import warn
import pdb

from mbuild.ff.opls_forcefield import OplsForceField
import mbuild.unit as units
from mbuild.rules import RuleEngine
from mbuild.prototype import Prototype
from mbuild.plot import Plot


class OplsRules(RuleEngine):
    """Loads a forcefield and parameterizes a compound."""
    def __init__(self, compound, force_field):
        super(OplsRules, self).__init__(compound)
        self.force_field = force_field

    def execute(self):
        unique_bond_types = set()

        for atom in self.compound.atoms():
            if atom.kind != "G":
                try:
                    unique_bond_types.add(Prototype.getAttr(atom.kind,
                            "bond_type"))
                    unique_bond_types.add(Prototype.getAttr(atom.kind,
                            "bond_type")[0]+"*")
                except:
                    warn('Found atomtype in component that is not present '
                            'in forcefield: {0}'.format(atom))

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

            # Find all atoms that match the bond type labels
            # (kind, alias with optional wildcards).
            atomTypes1 = self.force_field.findAtomTypes(pair[0])
            atomTypes2 = self.force_field.findAtomTypes(pair[1])


            # For every combination of matching atom kinds, create a bond type
            for (atomType1, atomType2) in product(atomTypes1, atomTypes2):
                if(atomType1 < atomType2):
                    pair = (atomType1, atomType2)
                else:
                    pair = (atomType1, atomType2)

                # Get all atoms that form this pair and try to add a bond.
                r_err = 0.2 * r 
                self.add_bond(pair[0], pair[1],
                        (r - 2*r_err)._value, (r + r_err)._value,
                        "{0}-{1}".format(pair[0], pair[1]))

if __name__ == "__main__":
    from mbuild.components.mpc_monomer import MpcMonomer
    m = MpcMonomer()

    ff = OplsForceField()
    ff.get_atom_types(m)

    rules = OplsRules(m, ff)
    rules.execute()

    Plot(m).show()


