from itertools import combinations, product
from warnings import warn
import pdb

from mbuild.ff.opls_forcefield import OplsForceField
import mbuild.unit as units
from mbuild.rules import RuleEngine
from mbuild.prototype import Prototype


class OplsRules(RuleEngine):
    """Loads a forcefield and parameterizes a compound."""

    def __init__(self, compound, force_field):
        super(OplsRules, self).__init__(compound)
        self.force_field = force_field

    def execute(self, verbose=False):
        pairs = [(i, j) for i,j in self.force_field.bond_types]

        for pair in pairs:
            r, k = self.force_field.bond_types[pair]
            if verbose:
                print "Found for pair {0}".format(pair)

            r_err = 0.2 * r
            self.add_bond(pair[0], pair[1],
                    (r - 2*r_err)._value, (r + r_err)._value,
                    "{0}-{1}".format(pair[0], pair[1]))

if __name__ == "__main__":
    from mbuild.components.mpc_monomer import MpcMonomer
    m = MpcMonomer()

    ff = OplsForceField().prune(m)

    # import pdb
    # pdb.set_trace()

    rules = OplsRules(m, ff)
    rules.execute()

    from mbuild.plot import Plot
    Plot(m).show()


