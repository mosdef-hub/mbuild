import os
import pdb

class OplsForceField(object):
    """A container class for the OPLS forcefield
    """

    def __init__(self):
        """Populate the database using files bundled with GROMACS
        """
        self.atom_types = dict()
        gmx_lib = os.getenv('GMXLIB', '')

        atom_type_file = os.path.join(gmx_lib, 'oplsaa.ff/ffnonbonded.itp')
        with open(atom_type_file, 'r') as f:
            for line in f:
                fields = line.split()
                if fields[0] in [';', '[']:
                    if fields[0][0] == '#':
                        next(f)
                self.atom_types[fields[0]] = fields[1:8]

        parsable_keywords = {'[ bondtypes ]': self.parse_bond_types}
        bond_type_file = os.path.join(gmx_lib, 'oplsaa.ff/ffbonded_processed.itp')
        self.bond_types = dict()
        with open(bond_type_file, 'r') as f:
            for line in f:
                if line.strip():
                    keyword = line.strip()
                    fields = line.split()
                    if fields[0] in [';']:
                        if fields[0][0] == '#':
                            next(f)
                    elif keyword in parsable_keywords:
                        parsable_keywords[keyword](f)

    def get_atom_types(self, compound):
        """
        """
        for atom in compound.atoms:
            params = self.atom_types[atom.kind]
            Prototype(atom.kind,
                    bond_type=params[0],
                    atomic_num=int(params[1]),
                    mass=float(params[2]),
                    charge=float(params[3]),
                    ptype=params[4],
                    sigma=float(params[5]),
                    epsilon=float(params[6]))

    def parse_bond_types(self, f):
        """
        """
        for line in f:
            if not line.strip():
                break
            fields = line.split()
            pair = (fields[0], fields[1])
            self.bond_types[pair] = [fields[3], fields[4]]

if __name__ == "__main__":
    ff = OplsForceField()
    pdb.set_trace()
