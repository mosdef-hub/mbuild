# important! make sure forcefield_unit_manager is loaded first!!!
import forcefield_unit_manager
from scimath.units import unit_manager
from scimath.units.unit_db import UnitDB
from scimath.units.SI import coulomb, second, kilogram
from scimath.units.mass import gram
from scimath.units.quantity import Quantity
from scimath.units.time import usec, msec

__author__ = 'sallai'

from traits.api import HasTraits, Str, Float, Int
from scimath.units.quantity_traits import QuantityTrait

print "Unit manager in AtomType: " +str(unit_manager)


class AtomType(HasTraits):
    kind = Str
    bondType = Str
    atomicNumber = Int()
    mass = QuantityTrait(0, kilogram, 'mass')
    # charge = QuantityTrait(0, coulomb, 'charge')
    time = QuantityTrait(0, msec, 'time')
    sigma = Float()
    epsilon = Float()

if __name__ == "__main__":

    # q = Quantity(0, units=usec)
    # print q.family_name
    # C = AtomType(kind='C')
    # print C.mass
    # print coulomb.label
    pass