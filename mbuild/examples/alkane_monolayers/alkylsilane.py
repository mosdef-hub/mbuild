__author__ = 'CTK'

from mbuild.coordinate_transform import equivalence_transform
from mbuild.compound import Compound

from mbuild.examples.alkane.alkane import Alkane
from mbuild.examples.alkane_monolayers.silane import Silane


class AlkylSilane(Compound):
    """
    """
    def __init__(self, chain_length):
        """
        """
        super(AlkylSilane, self).__init__()

        alkane = Alkane(chain_length, cap_end=False)
        self.add(alkane, 'alkane')
        silane = Silane()
        self.add(silane, 'silane')
        equivalence_transform(self.alkane, self.alkane.down, self.silane.up)

        # Hoist silane port to AlkylSilane level.
        self.add(silane.down, 'down', containment=False)

if __name__ == "__main__":
    alkyl_silane = AlkylSilane(5)

    from mbuild.plot import Plot
    Plot(alkyl_silane, bonds=True, verbose=False).show()

    #from mbuild.treeview import TreeView
    # TreeView(m).show()
