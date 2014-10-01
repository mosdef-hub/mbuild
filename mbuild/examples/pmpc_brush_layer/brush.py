from numpy import pi

from mbuild.compound import Compound
from mbuild.tools.polymer import Polymer
from mbuild.coordinate_transform import equivalence_transform
from mbuild.examples.ethane.methyl import Methyl
from mpc_monomer import MpcMonomer
from initiator import Initiator
from silane import Silane


class Brush(Compound):
    """ """
    def __init__(self, chain_length=4, alpha=pi/4):
        Compound.__init__(self)

        # Add parts
        self.add(Silane(), 'silane')
        self.add(Initiator(), 'initiator')
        self.add(Polymer(MpcMonomer(alpha=alpha), n=chain_length,
                         port_labels=("top_port", "bottom_port")), 'pmpc')
        self.add(Methyl(), 'methyl')

        equivalence_transform(self.initiator, self.initiator.bottom_port, self.silane.top_port)
        equivalence_transform(self.pmpc, self.pmpc.bottom_port, self.initiator.top_port)
        equivalence_transform(self.methyl, self.methyl.down, self.pmpc.top_port)

        # make self.port point to silane.bottom_port
        self.add(self.silane.bottom_port, label="port", containment=False)

if __name__ == "__main__":
    m = Brush(chain_length=20, alpha=pi/4)
    m = m.to_trajectory()

    m.save(filename='brush.xyz')

    #from mbuild.plot import Plot
    #Plot(m).show()
