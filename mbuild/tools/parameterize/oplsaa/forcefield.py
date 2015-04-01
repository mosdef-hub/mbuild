import os

import simtk.unit as units

from mbuild.tools.parameterize.forcefield import Forcefield, ForcefieldAtomtype

FILE_DIR = os.path.dirname(os.path.realpath(__file__))


class OPLSForcefield(Forcefield):
    """A container class for the OPLS forcefield."""

    def __init__(self, topology):
        """Populate the database using files bundled with GROMACS."""
        super(OPLSForcefield, self).__init__(topology)
        
        self.DEG = units.degrees
        self.RAD = units.radians

        # GROMACS specific units
        self.MASS = units.amu
        self.CHARGE = units.elementary_charge
        self.DIST = units.nanometers
        self.ENERGY = units.kilojoules_per_mole

        nonbonded_file = os.path.join(FILE_DIR, 'ffnonbonded_processed.itp')
        with open(nonbonded_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.split()
                    if fields[0][0] in [';', '[', '#']:
                        continue
                    self.add_atom_type(fields[0],
                                       fields[1],
                                       int(fields[2]),
                                       float(fields[3]) * self.MASS,
                                       float(fields[4]),
                                       #fields[5]  #  ignore ptype
                                       float(fields[6]) * self.DIST,
                                       float(fields[7]) * self.ENERGY)

        parsable_keywords = {'[ bondtypes ]': self.parse_bond_types,
                             '[ angletypes ]': self.parse_angle_types,
                             '[ dihedraltypes ]': self.parse_dihedral_types}
        bonded_file = os.path.join(FILE_DIR, 'ffbonded_processed.itp')
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
        """Read bond parameter information. """
        for line in f_bonded:
            if not line.strip():
                break
            fields = line.split()
            if fields[0][0] in [';', '#']:
                continue
            pair = tuple(sorted((fields[0], fields[1]), key=lambda x: x.lower()))
            if int(fields[2]) == 1:
                r = float(fields[3]) * self.DIST
                k = float(fields[4]) * self.ENERGY / (self.DIST * self.DIST)
            self.bondtypes[pair] = (r.in_units_of(units.angstroms), k)

    def parse_angle_types(self, f_bonded):
        """Read angle parameter information. """
        for line in f_bonded:
            if not line.strip():
                break
            fields = line.split()
            if fields[0][0] in [';', '#']:
                continue
            sorted_02 = sorted([fields[0], fields[2]], key=lambda x: x.lower())
            triplet = (sorted_02[0], fields[1], sorted_02[1])

            if int(fields[3]) == 1:
                theta = float(fields[4]) * self.DEG
                k = float(fields[5]) * self.ENERGY / (self.RAD * self.RAD)
            self.angletypes[triplet] = (theta, k)

    def parse_dihedral_types(self, f_bonded):
        """Read dihedral parameter information. """
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

            # Sort by inside pair
            if fields[0] == fields[3]:
                sorted_12 = sorted([fields[1], fields[2]], key=lambda x: x.lower())
                if fields[1] == sorted_12[0]:
                    quartet = (fields[0], fields[1], fields[2], fields[3])
                else:
                    quartet = (fields[3], fields[2], fields[1], fields[0])
            # Sort by outside pair
            else:
                sorted_03 = sorted([fields[0], fields[3]], key=lambda x: x.lower())
                if fields[0] == sorted_03[0]:
                    quartet = (fields[0], fields[1], fields[2], fields[3])
                else:
                    quartet = (fields[3], fields[2], fields[1], fields[0])
            self.dihedraltypes[quartet] = fields[5:end_of_params]

    def add_atom_type(self, opls_type, bond_type=None, atomic_number=0,
                      mass=0.0 * units.amu,
                      charge=0.0 * units.elementary_charge,
                      sigma=0.0 * units.angstroms,
                      epsilon=0.0 * units.kilojoules_per_mole):
        """
        """
        self.atomtypes[opls_type] = ForcefieldAtomtype(
            atomtype=opls_type, bondtype=bond_type, atomic_number=atomic_number,
            mass=mass, charge=charge, sigma=sigma, epsilon=epsilon)


