from mbuild.decorators import *
import mbuild.unit as units

class Atomtype(object):
    """
    """
    @accepts_compatible_units(None, None, None,
            units.amu,
            units.elementary_charge,
            units.angstroms,
            units.kilojoules_per_mole / units.angstroms**(2))
    def __init__(self, kind, alias=None, atomic_number=0,
            mass=0.0 * units.amu,
            charge=0.0 * units.elementary_charge,
            sigma=0.0 * units.angstroms,
            epsilon=0.0 * units.kilojoules_per_mole / units.angstroms**(2)):
        """
        """
        self.kind = kind
        self.alias = alias
        self.atomic_number = atomic_number
        self.mass = mass
        self.charge = charge
        self.sigma = sigma
        self.epsilon = epsilon

