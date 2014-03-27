import os
import pdb

from mbuild.prototype import Prototype
import mbuild.unit as units


class OplsForceField(object):
    """A container class for the OPLS forcefield."""

    def __init__(self):
        """Populate the database using files bundled with GROMACS."""
        self.atom_types = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.dihedral_types = dict()

        # GROMACS specific units
        self.DIST = units.nanometers
        self.ENERGY = units.kilojoules_per_mole

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
            if int(fields[2]) == 1:
                r = float(fields[3]) * self.DIST
                k = float(fields[3]) * self.ENERGY / (self.DIST * self.DIST)
            self.bond_types[pair] = (r.in_units_of(units.angstroms), k)

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

    def findAtomTypes(self, id):
        # if the id is the atom type, return the AtomType object
        if id in self.atom_types:
            return [self.atom_types[id]]

        matching_atom_types = []

        for at,bt in self.atom_types.iteritems():

            if str(id).endswith('*'):
                # id is a wildcard ending in *
                prefix = str(id)[:-1]

                if bt[0].startswith(prefix):
                    matching_atom_types.append(at)
                elif at.startswith(prefix):
                    matching_atom_types.append(at)
            else:
                # id is not a wildcard
                if id == bt[0]:
                    matching_atom_types.append(at)

        # return the matching_atom_types
        return matching_atom_types


if __name__ == "__main__":
    ff = OplsForceField()
    pdb.set_trace()
