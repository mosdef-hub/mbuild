from mbuild.coordinate_transform import equivalence_transform
from mbuild.atom import Atom
from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.tools.polymer import Polymer

from cco import Cco
from silane import Silane


class PegSilane(Compound):
    """
    """
    def __init__(self, n_peg):
        """
        """
        super(PegSilane, self).__init__()

        peg = Polymer(Cco(), n=n_peg, port_labels=('up', 'down'))
        self.add(peg, 'peg')

        hydrogen = Compound('H')
        hydrogen.add(Atom('H'), 'H')
        hydrogen.add(Port(anchor=hydrogen.H), 'down')
        self.add(hydrogen, 'hydrogen')
        equivalence_transform(self.peg, self.peg.up, self.hydrogen.down)

        silane = Silane()
        self.add(silane, 'silane')
        equivalence_transform(self.silane, self.silane.down, self.peg.down)

        # Hoist silane port to AlkylSilane level.
        self.add(silane.up, 'up', containment=False)

if __name__ == "__main__":
    pegsilane = PegSilane(10)
    pegsilane.visualize()

    #Plot(alkyl_silane, bonds=True, verbose=False).show()

    #from mbuild.treeview import TreeView
    # TreeView(m).show()
