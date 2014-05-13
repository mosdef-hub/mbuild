from itertools import combinations, product
from warnings import warn
import time
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
        """
        """
        start = time.time()
        pairs = self.force_field.bond_types.keys()
        for pair in pairs:
            r, k = self.force_field.bond_types[pair]
            if verbose:
                print "Found for pair {0}".format(pair)
            r_err = 0.2 * r
            self.add_bond(pair[0], pair[1],
                    (r - 2*r_err)._value, (r + r_err)._value,
                    "{0}-{1}".format(pair[0], pair[1]))
        print "    Done with bonds. ({0:.2f} s)".format(time.time() - start)

        start = time.time()
        triplets = self.force_field.angle_types.keys()
        for triplet in triplets:
            theta, k = self.force_field.angle_types[triplet]
            if verbose:
                print "Found for triplet {0}".format(triplet)
            self.add_angle(triplet[0], triplet[1], triplet[2],
                    "{0}-{1}-{2}".format(triplet[0], triplet[1], triplet[2]))
        print "    Done with angles. ({0:.2f} s)".format(time.time() - start)

        start = time.time()
        quartets = self.force_field.dihedral_types.keys()
        for quartet in quartets:
            cs = self.force_field.dihedral_types[quartet]
            if verbose:
                print "Found for quartet {0}".format(quartet)
            self.add_dihedral(quartet[0], quartet[1], quartet[2], quartet[3],
                    "{0}-{1}-{2}-{3}".format(quartet[0], quartet[1], quartet[2], quartet[3]))
        print "    Done with dihedrals. ({0:.2f} s)".format(time.time() - start)

if __name__ == "__main__":
    from mbuild.components.mpc_monomer import MpcMonomer
    m = MpcMonomer()

    ff = OplsForceField().prune(m)

    rules = OplsRules(m, ff)
    rules.execute()

    from mbuild.plot import Plot
    Plot(m).show()


