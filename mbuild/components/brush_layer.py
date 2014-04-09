from numpy import pi

from mbuild.plot import Plot
from mbuild.port import Port
from mbuild.compound import Compound

from mbuild.ff.opls_rules import OplsRules
from mbuild.ff.opls_forcefield import OplsForceField
import mbuild.unit as units

from mbuild.components.surface import Surface
from mbuild.components.initiator import Initiator
from mbuild.components.pmpc import Pmpc
from mbuild.components.silane import Silane
from mbuild.components.alkane_tail import AlkaneTail


class BrushLayer(Compound):
    """
    """

    def __init__(self, ctx={}, chain_length=4, alpha=pi/4, coverage=1):
        super(BrushLayer, self).__init__(ctx=ctx)
        self.add(Surface(ctx=ctx), 'surface')

        chains_on_surface = 0.0
        n_ports = sum(isinstance(part, Port) for part in self.surface.parts)
        for port in self.surface.parts:
            current_coverage = (chains_on_surface / n_ports ) * 100
            # Build a pMPC brush.
            if isinstance(port, Port) and current_coverage <  coverage:
                silane = Silane(ctx=ctx)
                silane.transform([(silane.bottom_port, port)])
                self.add(silane)

                initiator = Initiator(ctx=ctx)
                initiator.transform([(initiator.bottom_port, silane.top_port)])
                self.add(initiator)

                pmpc = Pmpc(ctx=ctx, n=chain_length, alpha=alpha)
                pmpc.transform([(pmpc.bottom_port, initiator.top_port)])
                self.add(pmpc)

                ch3 = AlkaneTail(ctx=ctx)
                ch3.transform([(ch3.female_port, pmpc.top_port)])
                self.add(ch3)

                chains_on_surface += 1

            elif current_coverage >= coverage:
                break

if __name__ == "__main__":
    m = BrushLayer(chain_length=1, alpha=pi/4, coverage=1)

    ff = OplsForceField()
    # TODO: add real parameters
    ff.add_atom_type(
            opls_type     = 'Si',
            bond_type     = 'Sisub',
            atomic_number = 14,
            mass          = 28.085 * units.amu,
            charge        = 0.84 * units.elementary_charge,
            sigma         = 3.5 * units.angstroms,
            epsilon       = 4.0 * units.kilojoules_per_mole)

    ff.add_atom_type(
            opls_type     = 'O',
            bond_type     = 'Osub',
            atomic_number = 8,
            mass          = 16.0 * units.amu,
            charge        = -0.42 * units.elementary_charge,
            sigma         = 3.5 * units.angstroms,
            epsilon       = 4.0 * units.kilojoules_per_mole)

    ff.add_atom_type(
            opls_type     = 'H',
            bond_type     = 'Hsub',
            atomic_number = 1,
            mass          = 1 * units.amu,
            charge        = 0.2 * units.elementary_charge,
            sigma         = 0.0 * units.angstroms,
            epsilon       = 0.0 * units.kilojoules_per_mole)

    ff.get_atom_types(m)
    rules = OplsRules(m, ff)
    rules.execute()

    Plot(m).show()
