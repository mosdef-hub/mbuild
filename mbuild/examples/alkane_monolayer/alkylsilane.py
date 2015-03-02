from __future__ import print_function

from mbuild.coordinate_transform import equivalence_transform
from mbuild.compound import Compound

from mbuild.examples.alkane.alkane import Alkane
from mbuild.components.small_groups.silane import Silane


class AlkylSilane(Compound):
    """A silane functionalized alkane chain with one Port. """
    def __init__(self, chain_length):
        super(AlkylSilane, self).__init__()

        alkane = Alkane(chain_length, cap_end=False)
        self.add(alkane, 'alkane')
        silane = Silane()
        self.add(silane, 'silane')
        equivalence_transform(self.alkane, self.alkane.down, self.silane.up)

        # Hoist silane port to AlkylSilane level.
        self.add(silane.down, 'down', containment=False)

if __name__ == "__main__":
    alkyl_silane = AlkylSilane(10)
    alkyl_silane = alkyl_silane.to_trajectory(show_ports=True)

    alkyl_silane.top.find_forcefield_terms()
    alkyl_silane.save('as.lammps', show_ports=True)
