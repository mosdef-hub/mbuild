import itertools
from mbuild.plot import Plot
from mbuild.prototype import Prototype
from mbuild.rules import RuleEngine
from wip.mpc_monolayer.mpc_monomer import MpcMonomer
# import wip.mpc_monolayer.mpc_monomer
from opls_forcefield import OplsForceField
import pdb
__author__ = 'sallai'


class OplsRules(RuleEngine):

    def __init__(self, compound, force_field):
        super(OplsRules, self).__init__(compound)
        self.force_field = force_field

    def execute(self):
        r_err = 0.3
        unique_bond_types = set()

        for atom in self.compound.atoms():
            if atom.kind != "G":
                unique_bond_types.add(Prototype.getAttr(atom.kind, "bond_type"))

        pairs = [ (i, j) for i,j  in itertools.combinations(unique_bond_types, 2)]


        for pair in pairs:

            if pair in self.force_field.bond_types:
                (r,k) = self.force_field.bond_types[pair]
            elif reversed(pair) in self.force_field.bond_types:
                (r,k) = self.force_field.bond_types[reversed(pair)]
            else:
                print "No bond type found for pair ("+pair[0] +","+pair[1]+")"
                continue

            # self.add_bond(atom1.bond_type, atom2.bond_type,
            #         r - r_err, r + r_err,
            #         "{0}-{1}".format(atom1.bond_type, atom2.bond_type))

if __name__ == "__main__":
    m = MpcMonomer()
    ff = OplsForceField()

    ff.get_atom_types(m)

    # Plot(m).show()

    r = OplsRules(m,ff)

    r.execute()



