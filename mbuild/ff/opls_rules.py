from itertools import combinations, product
from warnings import warn
import time
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

    def execute(self, verbose=False):
        unique_bond_types = set()

        now = time.time()
        for atom in self.compound.atoms():
            if atom.kind != "G":
                try:
                    unique_bond_types.add(Prototype.getAttr(atom.kind,
                            "bond_type"))
                    unique_bond_types.add("{0}*".format(
                            Prototype.getAttr(atom.kind, "bond_type")[0]))
                except:
                    warn('Found atomtype in component that is not present '
                            'in forcefield: {0}'.format(atom))

        pairs = [(i, j) for i,j in combinations(unique_bond_types, 2)]
        print "...Generated possible pairs: {0:.2f}".format(time.time() - now)

        now = time.time()
        for pair in pairs:
            if pair in self.force_field.bond_types:
                r, k = self.force_field.bond_types[pair]
                if verbose:
                    print "Found for pair {0}".format(pair)
            elif pair[::-1] in self.force_field.bond_types:
                r, k = self.force_field.bond_types[pair[::-1]]
                if verbose:
                    print "Found for pair {0}".format(pair)
            else:
                if verbose:
                    print "No bond type found for pair {0}".format(pair)
                continue

            # Find all atoms that match the bond type labels
            # (kind, alias with optional wildcards).
            atom_types_1 = self.force_field.find_atom_types(pair[0])
            pdb.set_trace()
            atom_types_2 = self.force_field.find_atom_types(pair[1])


            # For every combination of matching atom kinds, create a bond type
            for (atom_type_1, atom_type_2) in product(atom_types_1, atom_types_2):
                if atom_type_1 < atom_type_2:
                    pair = (atom_type_1, atom_type_2)
                else:
                    pair = (atom_type_1, atom_type_2)

                # Get all atoms that form this pair and try to add bonds.
                r_err = 0.2 * r 
                self.add_bond(pair[0], pair[1],
                        (r - 2*r_err)._value, (r + r_err)._value,
                        "{0}-{1}".format(pair[0], pair[1]))

        print "...Generated bonds {0:.2f}".format(time.time() - now)

if __name__ == "__main__":
    from mbuild.components.mpc_monomer import MpcMonomer
    m = MpcMonomer()

    ff = OplsForceField()
    ff.get_atom_types(m)

    rules = OplsRules(m, ff)
    rules.execute()

    Plot(m).show()


