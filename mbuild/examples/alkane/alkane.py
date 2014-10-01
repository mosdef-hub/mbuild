__author__ = 'sallai'

from mbuild.compound import Compound
from mbuild.tools.polymer import Polymer
from mbuild.coordinate_transform import equivalence_transform

from mbuild.examples.alkane.ch2 import Ch2
from mbuild.examples.ethane.methyl import Methyl


class Alkane(Compound):
    """ """
    def __init__(self, n=3, cap_front=True, cap_end=True):
        """Initialize an Alkane Compound.

        Args:
            n: Number of carbon atoms.
            cap_front: Add methyl group to beginning of chain ('down' port).
            cap_end: Add methyl group to end of chain ('up' port).
        """
        if n < 2:
            raise Exception('n must be 1 or more')
        Compound.__init__(self)

        # Adjust length of Polmyer for absence of methyl terminations.
        if not cap_front:
            n += 1
        if not cap_end:
            n += 1
        chain = Polymer(Ch2(), n=n-2, port_labels=('up', 'down'))
        self.add(chain, 'chain')

        if cap_front:
            self.add(Methyl(), "methyl_front")
            equivalence_transform(
                self.chain, self.chain.up, self.methyl_front.down)
        else:
            # Hoist port label to Alkane level.
            self.add(chain.up, 'up', containment=False)

        if cap_end:
            self.add(Methyl(), 'methyl_end')
            equivalence_transform(
                self.methyl_end, self.methyl_end.up, self.chain.down)
        else:
            self.add(chain.down, 'down', containment=False)


def main():
    n = 3
    alkane = Alkane(n=n, cap_front=True, cap_end=True)

    alkane = alkane.to_trajectory()
    alkane.top.find_forcefield_terms()


if __name__ == "__main__":
    main()
