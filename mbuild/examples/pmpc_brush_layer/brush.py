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
        super(Brush, self).__init__(self)

        # Add parts
        self.add(Silane(), 'silane')
        self.add(Initiator(), 'initiator')
        self.add(Polymer(MpcMonomer(alpha=alpha), n=chain_length,
                         port_labels=('up', 'down')), 'pmpc')
        self.add(Methyl(), 'methyl')

        equivalence_transform(self.initiator, self.initiator.down, self.silane.up)
        equivalence_transform(self.pmpc, self.pmpc.down, self.initiator.up)
        equivalence_transform(self.methyl, self.methyl.down, self.pmpc.up)

        # Make self.port point to silane.bottom_port
        self.add(self.silane.down, label="port", containment=False)

if __name__ == "__main__":
    m = Brush(chain_length=5, alpha=pi/4)
    m = m.to_trajectory()

    m.xyz[-1] *= 10.0
    m.save(filename='brush.mol2')

    #from mbuild.plot import Plot
    #Plot(m).show()
