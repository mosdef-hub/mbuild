import operator
from mbuild.plot import Plot
from mpc_monomer import MpcMonomer
from alkane_tail import AlkaneTail
from silane import Silane
from initiator import Initiator

__author__ = 'sallai'
from mbuild.compound import *
from mbuild.xyz import *
from mbuild.port import *

class Pmpc(Compound):


    def __init__(self, n=2, ctx={}):

        super(Pmpc, self).__init__(ctx=ctx)

        silane = Silane(ctx=ctx)
        self.add(silane, 'silane')
        initiator = Initiator(ctx=ctx)
        initiator.transform([(initiator.bottom_port, silane.top_port)])
        self.add(initiator, 'initiator')

        # n-1 times the body CH_2
        last_part = initiator
        for body_count in range(0, n):
            this_part = MpcMonomer(ctx=ctx, alpha=pi/3)
            self.add(this_part, 'monomer_'+str(body_count))
            if last_part is None:
                first_part = this_part
            else:
                this_part.transform([(this_part.bottom_port, last_part.top_port)])
            last_part = this_part

        ch3 = AlkaneTail(ctx=ctx)
        ch3.transform([(ch3.bottom_port, last_part.top_port)])
        self.add(ch3, "CH3")

        #top_port = Port()
        #translateTo = map(operator.sub,last_part.top_port.middle.pos, top_port.middle.pos)
        #top_port.transform(Translation(tuple(translateTo)))
        #self.add(top_port, "top_port")

        bottom_port = Port()
        #translateTo = map(operator.sub, silane.bottom_port.middle.pos, bottom_port.middle.pos)
        #bottom_port.transform(Translation(tuple(translateTo)))
        bottom_port.transform([(bottom_port, silane.bottom_port)])
        self.add(bottom_port, "bottom_port")

if __name__ == "__main__":
    m = Pmpc(n=50)
    # TreeView(m, verbose=True).show()
    Plot(m, verbose=True).show()
