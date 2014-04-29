from numpy import pi

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


class Brush(Compound):
    """
    """

    def __init__(self, ctx={}, chain_length=4, alpha=pi/4):
        super(Brush, self).__init__(ctx=ctx)

        silane = Silane(ctx=ctx)
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

        # make self.port point to silane.bottom_port
        self.add(silane.bottom_port, label="port", containment=False)


if __name__ == "__main__":
    m = Brush(chain_length=5, alpha=pi/4)

    from mbuild.plot import Plot

    Plot(m).show()
