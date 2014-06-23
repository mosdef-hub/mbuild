import mbuild.unit as units
from examples.ethane.ethane import Ethane
from mbuild.compound import Compound
from mbuild.ff.opls_forcefield import OplsForceField
from mbuild.ff.opls_rules import OplsRules
from mbuild.plot import Plot
from mbuild.prototype import Prototype
from mbuild.unit.quantity import Quantity
from mbuild.xyz import Xyz
import os
import time
import sys
__author__ = 'sallai'


class WaterBox(Compound):
    """
    """

    def __init__(self, ctx=None):
        """
        """
        if not ctx:
            ctx = {}

        super(WaterBox, self).__init__(ctx=ctx)

        oplsff = OplsForceField()

        oplsff.add_atom_type(
            opls_type     = 'OW',
            bond_type     = 'OW',
            atomic_number = 8,
            mass          = 16.0 * units.amu,
            charge        = -0.42 * units.elementary_charge,
            sigma         = 3.5 * units.angstroms,
            epsilon       = 4.0 * units.kilojoules_per_mole)

        oplsff.add_atom_type(
            opls_type     = 'HW1',
            bond_type     = 'HW',
            atomic_number = 1,
            mass          = 1.0 * units.amu,
            charge        = 0.0  * units.elementary_charge,
            sigma         = 0 * units.angstroms,
            epsilon       = 1.0 * units.kilojoules_per_mole)

        oplsff.add_atom_type(
            opls_type     = 'HW2',
            bond_type     = 'HW',
            atomic_number = 1,
            mass          = 1.0 * units.amu,
            charge        = 0.0  * units.elementary_charge,
            sigma         = 0 * units.angstroms,
            epsilon       = 1.0 * units.kilojoules_per_mole)

        oplsff.add_atom_type(
            opls_type     = 'H',
            bond_type     = 'HC',
            atomic_number = 1,
            mass          = 1.0 * units.amu,
            charge        = 0.0  * units.elementary_charge,
            sigma         = 0 * units.angstroms,
            epsilon       = 1.0 * units.kilojoules_per_mole)


        oplsff.add_atom_type(
            opls_type     = 'C',
            bond_type     = 'CT',
            atomic_number = 6,
            mass          = 12.0 * units.amu,
            charge        = 0.0  * units.elementary_charge,
            sigma         = 0 * units.angstroms,
            epsilon       = 1.0 * units.kilojoules_per_mole)



        # load water
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        xyz_path = os.path.join(current_dir, 'spc216.xyz')

        w = Xyz(xyz_path, labels=False)
        self.add(w, 'water')
        wff = oplsff.prune(w)
        wff.init_prototypes()
        print "Generating water topology..."
        start = time.time()
        rules = OplsRules(w, wff)
        rules.execute(verbose=False)
        print "Done. ({0:.2f} s)".format(time.time() - start)

        # create ethane
        ethane = Ethane()
        self.add(ethane, "ethane")
        ethaneff = oplsff.prune(ethane)
        ethaneff.init_prototypes()
        print "Generating ethane topology..."
        start = time.time()
        rules = OplsRules(w, ethaneff)
        rules.execute(verbose=False)
        print "Done. ({0:.2f} s)".format(time.time() - start)

        # Prototype('OW', sigma=3.15)
        # Prototype('C', sigma=3.5)
        # Prototype('H', sigma=2.5)


        for o in self.getAtomListByKind(kind='OW'):
            o_sigma = Prototype.getAttr(o.kind, "sigma")
            if isinstance(o_sigma, Quantity):
                o_sigma = o_sigma._value
            print "o_sigma=" + str(o_sigma)
            for neighbor in ethane.getAtomsInRange(o.pos, 10, maxItems=10, kind='*'):
                if neighbor is o:
                    continue
                neighbor_sigma = Prototype.getAttr(neighbor.kind, "sigma", default=float("-inf"))
                if isinstance(neighbor_sigma, Quantity):
                    neighbor_sigma = neighbor_sigma._value
                print "neighbor_sigma=" + str(neighbor_sigma)
                if o.distance(neighbor) < o_sigma + neighbor_sigma:
                    print str(o) + " is close to " + str(neighbor)
                    print "removing " + str(o)

                    for b in o.bonds:
                        print str(b)
                        self.remove(b.atom1)
                        self.remove(b.atom2)
                        print "removing " + str(b.atom1)
                        print "removing " + str(b.atom2)

if __name__ == "__main__":
    wb = WaterBox()

    # for a1 in wb.atoms():
    #     for a2 in wb.atoms():
    #         d = a1.distance(a2)
    #         if d < 10:
    #             print d

    Plot(wb, bonds=True).show()