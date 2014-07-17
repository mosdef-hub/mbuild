from numpy import pi
from mbuild.coordinate_transform import Translation, transform

from mbuild.port import Port
from mbuild.compound import Compound

from initiator import Initiator
from pmpc import Pmpc
from silane import Silane
from methyl import Methyl


class Brush(Compound):
    """
    """

    def __init__(self, chain_length=4, alpha=pi/4):
        Compound.__init__(self)

        silane = Silane()
        self.add(silane)

        initiator = Initiator()
        transform(initiator, [(initiator.labels['bottom_port'], silane.labels['top_port'])])
        self.add(initiator)

        pmpc = Pmpc(n=chain_length, alpha=alpha)
        transform(pmpc, [(pmpc.labels['bottom_port'], initiator.labels['top_port'])])
        self.add(pmpc)

        ch3 = Methyl()
        transform([ch3, (ch3.labels['female_port'], pmpc.labels['top_port'])])
        self.add(ch3)

        # make self.port point to silane.bottom_port
        self.add(silane.labels['bottom_port'], label="port", containment=False)


if __name__ == "__main__":
    m = Brush(chain_length=5, alpha=pi/4)

    from mbuild.plot import Plot

    Plot(m).show()
