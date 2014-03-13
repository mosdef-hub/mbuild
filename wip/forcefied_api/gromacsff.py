from scimath.units.mass import gram
from scimath.units.quantity import Quantity
from AtomType import AtomType

__author__ = 'sallai'
import os

from traits.api import HasTraits, Str, Float
from scimath.units.quantity_traits import QuantityTrait



class GromacsFF(object):
    """A reader for GROMACS forcefield files
    """

    def __init__(self, forcefield_dir):

        self.atom_types = dict()

        gmx_lib = os.getenv('GMXLIB', 'gromacs_forcefields')

        atom_type_file = os.path.join(gmx_lib, forcefield_dir, 'ffnonbonded.itp')
        with open(atom_type_file, 'r') as f:
            for line in f:
                fields = line.split()
                if not fields[0][0] in [';', '[', '#']:
                    atomType = AtomType(kind=fields[0], bondType=fields[1], atomicNumber=fields[2], mass=Quantity(fields[3], units=gram, family_name='mass') , charge=fields[4], sigma=fields[6], epsilon=fields[7])
                    # self.atom_types[fields[0]] = fields[1:8]
                    print atomType

        parsable_keywords = {'[ bondtypes ]': self.parse_bond_types}
        bond_type_file = os.path.join(gmx_lib, forcefield_dir, 'ffbonded.itp')
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
    ff = GromacsFF("oplsaa.ff")

