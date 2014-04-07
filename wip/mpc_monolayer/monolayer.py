
__author__ = 'sallai'

import pdb

from mbuild.plot import Plot
from mbuild.port import Port
from mbuild.treeview import TreeView
from mbuild.compound import *
#from wip.forcefield_rules.opls_rules import OplsRules
#from wip.forcefield_rules.opls_forcefield import OplsForceField
from mbuild.ff.opls_rules import OplsRules
from mbuild.ff.opls_forcefield import OplsForceField
import mbuild.unit as units

from surface import Surface
from initiator import Initiator
from pmpc import Pmpc
from silane import Silane
from alkane_tail import AlkaneTail

class Monolayer(Compound):

    def __init__(self, chain_length=4, coverage=13, ctx={}):

        super(Monolayer, self).__init__(ctx=ctx)

        self.add(Surface(ctx=ctx),'surface')

        coverage_cnt = 0
        for port in self.surface.parts:
            if isinstance(port, Port) and (coverage_cnt % coverage) == 0:
                print coverage_cnt
                silane = Silane(ctx=ctx)
                silane.transform([(silane.bottom_port, port)])
                self.add(silane)

                initiator = Initiator(ctx=ctx)
                initiator.transform([(initiator.bottom_port, silane.top_port)])
                self.add(initiator)

                pmpc = Pmpc(ctx=ctx, n=chain_length)
                pmpc.transform([(pmpc.bottom_port, initiator.top_port)])
                self.add(pmpc)

                ch3 = AlkaneTail(ctx=ctx)
                ch3.transform([(ch3.female_port, pmpc.top_port)])
                self.add(ch3)

            coverage_cnt = coverage_cnt+1

if __name__ == "__main__":
    m = Monolayer(chain_length=13, coverage=100)
    ff = OplsForceField()

    ff.add_atom_type(opls_type='Si',
            bond_type='Sisub',
            atomic_number=14,
            mass=28.085 * units.amu,
            charge=0.84 * units.elementary_charge,
            sigma=3.5 * units.angstroms,
            epsilon=4.0 * units.kilojoules_per_mole / units.nanometers**(2))

    ff.add_atom_type(opls_type='O',
            bond_type='Osub',
            atomic_number=8,
            mass=16.0 * units.amu,
            charge=-0.42 * units.elementary_charge,
            sigma=3.5 * units.angstroms,
            epsilon=4.0 * units.kilojoules_per_mole / units.nanometers**(2))

    ff.add_atom_type(opls_type='H',
            bond_type='Hsub',
            atomic_number=1,
            mass=1 * units.amu,
            charge=0.2 * units.elementary_charge,
            sigma=0.0 * units.angstroms,
            epsilon=0.0 * units.kilojoules_per_mole / units.nanometers**(2))


    ff.get_atom_types(m)
    rules = OplsRules(m, ff)
    rules.execute()
    Plot(m).show()
    # print [(label,atom.pos) for label, atom in m.atoms()]
    #TreeView(m).show()
    # m.plot(labels=False, verbose=True)
