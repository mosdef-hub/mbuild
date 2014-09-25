__author__ = 'sallai'

from mbuild.compound import Compound
from mbuild.examples.ethane.methyl import Methyl
from mbuild.coordinate_transform import equivalence_transform


class Ethane(Compound):
    """ """

    def __init__(self):
        super(Ethane, self).__init__(kind='Ethane')
        self.add(Methyl(), "m1")
        self.add(Methyl(), "m2")

        equivalence_transform(self.m1, self.m1.up, self.m2.down)

if __name__ == "__main__":
    ethane = Ethane()
    
    ethane.save('ethane.pdb')
    ethane.save('ethane.hoomdxml')
    
    
    ethane = ethane.to_trajectory()
    ethane.top.find_forcefield_terms()

    print len(ethane.top._ff_bonds)
    print len(ethane.top._ff_angles)
    print len(ethane.top._ff_dihedrals)

    #from mbuild.plot import Plot
    #Plot(ethane, verbose=False, atoms=True, bonds=True, angles=False, dihedrals=False).show()

    # from mbuild.treeview import TreeView
    # TreeView(ethane).show()
