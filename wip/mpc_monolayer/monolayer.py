
__author__ = 'sallai'

from mbuild.port import Port
from mbuild.treeview import TreeView
from mbuild.compound import *
from surface import Surface
from initiator import Initiator
from pmpc import Pmpc
from silane import Silane
from alkane_tail import AlkaneTail

class Monolayer(Compound):

    def __init__(self, chain_length=4, coverage=13, ctx={}):

        super(Monolayer, self).__init__(ctx=ctx)

        self.add(Surface(ctx=ctx),'surface')

        coverage_cnt = 0
        for port in self.surface.parts:
            if isinstance(port, Port) and (coverage_cnt % coverage) == 0:
                print coverage_cnt
                silane = Silane(ctx=ctx)
                silane.transform([(silane.bottom_port, port)])
                self.add(silane)

                initiator = Initiator(ctx=ctx)
                initiator.transform([(initiator.bottom_port, silane.top_port)])
                self.add(initiator)

                pmpc = Pmpc(ctx=ctx, n=chain_length)
                pmpc.transform([(pmpc.bottom_port, initiator.top_port)])
                self.add(pmpc)

                ch3 = AlkaneTail(ctx=ctx)
                ch3.transform([(ch3.female_port, pmpc.top_port)])
                self.add(ch3)

            coverage_cnt = coverage_cnt+1

if __name__ == "__main__":
    m = Monolayer(chain_length=13)
    # print [(label,atom.pos) for label, atom in m.atoms()]
    TreeView(m).show()
    # m.plot(labels=False, verbose=True)