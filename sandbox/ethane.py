from mbuild.moleculemodel import MoleculeModel
from mbuild.mpcchainrules import MpcChainRules
from mbuild.plot import Plot
from mbuild.rules import RuleEngine
from mbuild.treeview import TreeView

__author__ = 'sallai'

from alkane_tail import AlkaneTail
from mbuild.compound import *


class Ethane(Compound):

    def __init__(self, ctx={}):
        super(Ethane, self).__init__(kind='Ethane', ctx=ctx)
        # two tails
        self.add(AlkaneTail(), 'top_tail')
        self.add(AlkaneTail(), 'bottom_tail')
        # transform bottom_tail to top_tail's coordinate system with point equivalencies
        # m.bottom_tail.transform(
        #                     [
        #                             (m.bottom_tail.female_port.top, m.top_tail.male_port.top),
        #                             (m.bottom_tail.female_port.left, m.top_tail.male_port.left),
        #                             (m.bottom_tail.female_port.right, m.top_tail.male_port.right)
        #                     ])

        # transform bottom_tail to top_tail's coordinate system with subcomponent equivalencies
        self.bottom_tail.transform([(self.bottom_tail.female_port, self.top_tail.male_port)])


# class MyRules(RuleEngine):
#
#     @classmethod
#     def create(cls, model):
#         re = super(MyRules, cls).create(model)
#         return re
#
#
#     def execute(self):
#         self.add_bond(C, H, .9, 1.1, "c-h", (1, 1, 1))
#         self.add_bond(C, C, .9, 2.1, "c-c", (0, 1, 0))
#
#         self.add_angle(H, C, C, "h-c-c", color=(1, 0.5, 0))
#
#         self.add_dihedral(H, C, C, H, "h-c-c-h", color=(0, 0, 0))
#
#


if __name__ == "__main__":
    ethane = Ethane()
    print ethane.bottom_tail.c

    # Plot(ethane).show()
    TreeView(ethane).show()
