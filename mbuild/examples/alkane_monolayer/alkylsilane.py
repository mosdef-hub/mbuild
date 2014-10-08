__author__ = 'CTK'

from mbuild.coordinate_transform import equivalence_transform
from mbuild.compound import Compound

from mbuild.examples.alkane.alkane import Alkane
from silane import Silane


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
    alkyl_silane = AlkylSilane(100)
    alkyl_silane = alkyl_silane.to_trajectory()

    import cProfile, pstats, StringIO
    pr = cProfile.Profile()
    pr.enable()

    alkyl_silane.top.find_forcefield_terms()

    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()

    print alkyl_silane.n_atoms
    print alkyl_silane.top.n_ff_bonds
    print alkyl_silane.top.n_ff_angles
    print alkyl_silane.top.n_ff_dihedrals

    #from mbuild.plot import Plot
    #Plot(alkyl_silane, bonds=True, verbose=False).show()

    #from mbuild.treeview import TreeView
    # TreeView(m).show()
