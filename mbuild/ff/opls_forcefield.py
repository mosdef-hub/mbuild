from itertools import product
import os
import pdb

import numpy as np

import mbuild
from mbuild.ff.forcefield import ForceField
from mbuild.prototype import Prototype
import mbuild.unit as units
from mbuild.decorators import *
from mbuild.ff.atom_type import Atomtype


class OplsForceField(ForceField):
    """A container class for the OPLS forcefield."""

    def __init__(self):
        """Populate the database using files bundled with GROMACS."""
        super(OplsForceField, self).__init__(name='opls')

        self.DEG = units.degrees
        self.RAD = units.radians

        # GROMACS specific units
        self.MASS = units.amu
        self.CHARGE = units.elementary_charge
        self.DIST = units.nanometers
        self.ENERGY = units.kilojoules_per_mole

        nonbonded_file = os.path.join(mbuild.__path__[0], 'ff',
                'ffnonbonded_processed.itp')
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
                            float(fields[4]) * self.CHARGE,
                            #fields[5]  #  ignore ptype
                            float(fields[6]) * self.DIST,
                            float(fields[7]) * self.ENERGY)

        parsable_keywords = {'[ bondtypes ]': self.parse_bond_types,
                '[ angletypes ]': self.parse_angle_types,
                '[ dihedraltypes ]': self.parse_dihedral_types}
        bonded_file = os.path.join(mbuild.__path__[0], 'ff',
                'ffbonded_processed.itp')
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
                r = r.in_units_of(units.angstroms)
                k = float(fields[4]) * self.ENERGY / (self.DIST * self.DIST)
                self.add_bond_type(pair, r, k)
            else:
                raise Exception("Found non-opls angle in forcefield")

    def parse_angle_types(self, f_bonded):
        """Read angle parameter information."""
        for line in f_bonded:
            if not line.strip():
                break
            fields = line.split()
            if fields[0][0] in [';', '#']:
                continue
            triplet = (fields[0], fields[1], fields[2])
            if int(fields[3]) == 1:
                theta = float(fields[4]) * self.DEG
                k = float(fields[5]) * self.ENERGY / (self.RAD * self.RAD)
                self.add_angle_type(triplet, theta, k)
            else:
                raise Exception("Found non-opls angle in forcefield")

    def parse_dihedral_types(self, f_bonded):
        """Read dihedral parameter information."""
        for line in f_bonded:
            if not line.strip():
                break
            fields = line.split()
            if fields[0][0] in [';', '#']:
                continue
            quartet = (fields[0], fields[1], fields[2], fields[3])
            # TODO: convert all incoming dihedrals to OPLS style
            if fields[4].isdigit():
                if int(fields[4]) == 3:
                    cs = map(float, fields[5:11])
                    for i, c in enumerate(cs):
                        cs[i] = c * self.ENERGY
                    self.add_dihedral_type(quartet, *cs)
                elif int(fields[4]) == 4:
                    print "Found improper. Ignoring"
                    #params = np.array(fields[5:11], dtype=float)
                    #self.dihedral_types[quartet] = fields[5:end_of_params]

    @accepts_compatible_units(None, None, None,
            units.amu,
            units.elementary_charge,
            units.angstroms,
            units.kilojoules_per_mole)
    def add_atom_type(self, opls_type, bond_type=None, atomic_number=0,
            mass=0.0 * units.amu,
            charge=0.0 * units.elementary_charge,
            sigma=0.0 * units.angstroms,
            epsilon=0.0 * units.kilojoules_per_mole):
        """
        """
        self.atom_types[opls_type] = Atomtype(kind=opls_type, alias=bond_type,
                atomic_number=atomic_number, mass=mass, charge=charge,
                sigma=sigma, epsilon=epsilon)

    @accepts_compatible_units(None,
            units.angstroms,
            units.kilojoules_per_mole * units.nanometers**(-2))
    def add_bond_type(self, pair, r, k):
        self.bond_types[pair] = (r, k)

    @accepts_compatible_units(None,
            units.degrees,
            units.kilojoules_per_mole * units.radians**(-2))
    def add_angle_type(self, triplet, theta, k):
        self.angle_types[triplet] = (theta, k)

    @accepts_compatible_units(None,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole)
    def add_dihedral_type(self, quartet, c0, c1, c2, c3, c4, c5):
        self.dihedral_types[quartet] = (c0, c1, c2, c3, c4, c5)

if __name__ == "__main__":
    ff = OplsForceField()
    pdb.set_trace()
