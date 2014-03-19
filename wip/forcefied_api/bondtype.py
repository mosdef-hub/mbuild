# important! make sure forcefield_unit_manager is loaded first!!!
from traits.trait_types import Tuple
import forcefield_unit_manager
from scimath.units import unit_manager
from scimath.units.unit_db import UnitDB
from scimath.units.SI import coulomb, second, kilogram
from scimath.units.mass import gram
from scimath.units.quantity import Quantity
from scimath.units.time import usec, msec

__author__ = 'sallai'

from traits.api import HasTraits, Str, CFloat, Int
from scimath.units.quantity_traits import QuantityTrait


class BondType(HasTraits):
    bond_types = Tuple(('',''))
    r = CFloat('nan')
    k = CFloat('nan')


if __name__ == "__main__":
    pass