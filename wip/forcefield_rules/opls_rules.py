from mbuild.rules import RuleEngine
from mpc_monomer import MpcMonomer

__author__ = 'sallai'

class OplsRules(RuleEngine):

    def execute(self):
        r_err = 0.3
        for atom in atoms:
            unique_bond_types.add(atom.bond_type)

        pairs = filter_combo(unique_bond_types)
        for atom1, atom2 in pairs:

            r, k = fetch_bondtype(atom1.bond_type, atom2.bond_type)
            self.add_bond(atom1.bond_type, atom2.bond_type,
                    r - r_err, r + r_err, 
                    "{0}-{1}".format(atom1.bond_type, atom2.bond_type)

if __name__ == "__main__":
    m = MpcMonomer()
    m.pull_atomtypes("opls")
    r = OplsRules(m)
    r.execute()



