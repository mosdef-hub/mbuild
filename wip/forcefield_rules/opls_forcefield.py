import os
import pdb
from mbuild.prototype import Prototype


class OplsForceField(object):
    """A container class for the OPLS forcefield
    """

    def __init__(self):
        """Populate the database using files bundled with GROMACS
        """
        self.atom_types = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.dihedral_types = dict()

        #gmx_lib = os.getenv('GMXLIB', 'gromacs_forcefields')
        #nonbonded_file = os.path.join("../forcefied_api/gromacs_forcefields",
        #       'oplsaa.ff/ffnonbonded.itp')
        nonbonded_file = "ffnonbonded_processed.itp"
        with open(nonbonded_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.split()
                    if fields[0][0] in [';', '[', '#']:
                        continue
                    self.atom_types[fields[0]] = fields[1:8]

        parsable_keywords = {'[ bondtypes ]': self.parse_bond_types,
                '[ angletypes ]': self.parse_angle_types,
                '[ dihedraltypes ]': self.parse_dihedral_types}
        #bonded_file = os.path.join("../forcefied_api/gromacs_forcefields",
        #        'oplsaa.ff/ffbonded.itp')
        bonded_file = "ffbonded_processed.itp"
        with open(bonded_file, 'r') as f_bonded:
            for line in f_bonded:
                if line.strip():
                    keyword = line.strip()
                    fields = line.split()
                    if fields[0][0] in [';', '#']:
                        continue
                    elif keyword in parsable_keywords:
                        parsable_keywords[keyword](f_bonded)

    def parse_bond_types(self, f_bonded):
        """Read bond parameter information."""
        for line in f_bonded:
            if not line.strip():
                break
            fields = line.split()
            if fields[0][0] in [';', '#']:
                continue
            pair = (fields[0], fields[1])
            self.bond_types[pair] = [float(field) for field in fields[3:5]]

    def parse_angle_types(self, f_bonded):
        """Read angle parameter information."""
        for line in f_bonded:
            if not line.strip():
                break
            fields = line.split()
            if fields[0][0] in [';', '#']:
                continue
            triplet = (fields[0], fields[1], fields[2])
            self.angle_types[triplet] = [float(field) for field in fields[4:6]]

    def parse_dihedral_types(self, f_bonded):
        """Read dihedral parameter information."""
        for line in f_bonded:
            if not line.strip():
                break
            fields = line.split()
            if fields[0][0] in [';', '#']:
                continue
            if ';' in fields:
                end_of_params = fields.index(';')
            else:
                end_of_params = len(fields)

            quartet = (fields[0], fields[1], fields[2], fields[3])
            self.dihedral_types[quartet] = fields[5:end_of_params]



    def get_atom_types(self, compound):
        """Parameterize a compound with atom specific information."""
        for atom in compound.atoms():
            if atom.kind in self.atom_types:
                params = self.atom_types[atom.kind]
                Prototype(atom.kind,
                        bond_type=params[0],
                        atomic_num=int(params[1]),
                        mass=float(params[2]),
                        charge=float(params[3]),
                        ptype=params[4],
                        sigma=float(params[5]),
                        epsilon=float(params[6]))

if __name__ == "__main__":
    ff = OplsForceField()
    pdb.set_trace()
