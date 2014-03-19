# important! make sure forcefield_unit_manager is loaded first!!!
import forcefield_unit_manager
from scimath.units import unit_manager
from scimath.units.unit_db import UnitDB
from scimath.units.SI import coulomb, second, kilogram
from scimath.units.mass import gram
from scimath.units.quantity import Quantity
from scimath.units.time import usec, msec

__author__ = 'sallai'

from traits.api import HasTraits, Str, CFloat, CInt
from scimath.units.quantity_traits import QuantityTrait


class AtomType(HasTraits):
    kind = Str
    alias = Str
    atomicNumber = CInt(0)
    mass = QuantityTrait(float('nan'), kilogram, 'mass')
    charge = QuantityTrait(float('nan'), coulomb, 'charge')
    sigma = CFloat('nan')
    epsilon = CFloat('nan')

if __name__ == "__main__":
    pass