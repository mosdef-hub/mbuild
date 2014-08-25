__author__ = 'sallai'

import operator

from mpc_monomer import MpcMonomer
from mbuild.file_formats.xyz import *
from mbuild.port import *


class Pmpc(Compound):


    def __init__(self, n=2, ctx={}):

        if n < 1:
            raise Exception('n must be 1 or more')

        super(Pmpc, self).__init__(ctx=ctx)


        # n-1 times the body CH_2
        last_part = None
        for body_count in range(0, n):
            this_part = MpcMonomer(ctx=ctx, alpha=pi/6)
            self.add(this_part, 'monomer_'+str(body_count))
            if last_part is None:
                first_part = this_part
            else:
                this_part.transform([(this_part.bottom_port, last_part.top_port)])
            last_part = this_part

        top_port = Port()
        top_port.transform(RotationAroundZ(pi))
        translateTo = map(operator.sub,last_part.top_port.middle.pos, top_port.middle.pos)
        top_port.transform(Translation(tuple(translateTo)))
        self.add(top_port, "top_port")

        bottom_port = Port()
        bottom_port.transform(RotationAroundZ(pi))
        translateTo = map(operator.sub,first_part.bottom_port.middle.pos, bottom_port.middle.pos)
        bottom_port.transform(Translation(tuple(translateTo)))
        self.add(bottom_port, "bottom_port")



if __name__ == "__main__":
    m = Pmpc(n=13)
    # TreeView(m, verbose=True).show()
    Plot(m).show()
