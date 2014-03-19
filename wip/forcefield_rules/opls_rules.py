__author__ = 'sallai'

from itertools import combinations
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
        r_err = 0.3 * units.angstroms
        unique_bond_types = set()

        for atom in self.compound.atoms():
            if atom.kind != "G":
                unique_bond_types.add(Prototype.getAttr(atom.kind,
                        "bond_type"))

        pairs = [(i, j) for i,j in combinations(unique_bond_types, 2)]

        for pair in pairs:
            if pair in self.force_field.bond_types:
                r, k = self.force_field.bond_types[pair]
            elif pair[::-1] in self.force_field.bond_types:
                r, k = self.force_field.bond_types[pair[::-1]]
            else:
                print "No bond type found for pair {0}".format(pair)
                continue

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
    pdb.set_trace()
    Plot(m).show()


