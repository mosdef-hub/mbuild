__author__ = 'sallai'

# from mbuild.plot import Plot
# from mbuild.rules import RuleEngine
from mbuild.treeview import TreeView
from alkane_tail import AlkaneTail
from mbuild.compound import *


class Ethane(Compound):

    def __init__(self, ctx={}):
        super(Ethane, self).__init__(kind='Ethane', ctx=ctx)
        # two tails
        self.add(AlkaneTail(), 'top_tail')
        self.add(AlkaneTail(), 'bottom_tail')

        # transform bottom_tail to top_tail's coordinate system with point equivalencies
        self.bottom_tail.transform(
                            [
                                    (self.bottom_tail.female_port.top, self.top_tail.male_port.top),
                                    (self.bottom_tail.female_port.middle, self.top_tail.male_port.middle),
                                    (self.bottom_tail.female_port.left, self.top_tail.male_port.left),
                                    (self.bottom_tail.female_port.right, self.top_tail.male_port.right)
                            ])

        # # transform bottom_tail to top_tail's coordinate system with subcomponent equivalencies
        # self.bottom_tail.transform([(self.bottom_tail.female_port, self.top_tail.male_port)])

if __name__ == "__main__":
    ethane = Ethane()

    # Plot(ethane).show()
    TreeView(ethane).show()
