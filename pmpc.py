import operator
from mpc_monomer import MpcMonomer

__author__ = 'sallai'
from mbuild.compound import *
from mbuild.xyz import *
from mbuild.port import *

class Pmpc(Compound):


    def __init__(self, n=2, ctx={}):

        if n < 1:
            raise Exception('n must be 1 or more')

        super(Pmpc, self).__init__(ctx=ctx)


        # n-1 times the body CH_2
        last_part = None
        for body_count in range(0, n):
            this_part = MpcMonomer(ctx=ctx)
            self.add(this_part, 'monomer_'+str(body_count))
            if not last_part is None:
                this_part.transform([(this_part.bottom_port, last_part.top_port)])
            last_part = this_part

        top_port = Port()
        translateTo = map(operator.sub,last_part.top_port.top.pos, top_port.zero.pos)
        top_port.transform(Translation(tuple(translateTo)))
        self.add(top_port, "top_port")

        bottom_port = Port()
        translateTo = map(operator.sub,last_part.top_port.top.pos, bottom_port.zero.pos)
        bottom_port.transform(Translation(tuple(translateTo)))
        self.add(bottom_port, "bottom_port")



if __name__ == "__main__":
    m = Pmpc()
    TreeView(m, verbose=True).show()
