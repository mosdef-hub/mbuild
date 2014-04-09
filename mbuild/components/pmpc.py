import operator
import pdb

from numpy import pi

from mbuild.coordinate_transform import *
from mbuild.xyz import Xyz
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.plot import Plot

from mbuild.components.mpc_monomer import MpcMonomer


class Pmpc(Compound):

    def __init__(self, ctx={}, n=2, alpha=pi/6):
        if n < 1:
            raise Exception('n must be 1 or more')
        super(Pmpc, self).__init__(ctx=ctx)

        # n-1 times the body CH_2
        last_part = None
        for body_count in range(0, n):
            this_part = MpcMonomer(ctx=ctx, alpha=alpha)
            self.add(this_part, 'monomer_{0}'.format(body_count))
            if last_part is None:
                first_part = this_part
            else:
                this_part.transform([(this_part.bottom_port,
                                      last_part.top_port)])
            last_part = this_part

        top_port = Port()
        top_port.transform(RotationAroundZ(pi))
        translate_to = map(operator.sub,
            last_part.top_port.middle.pos, top_port.middle.pos)
        top_port.transform(Translation(translate_to))
        self.add(top_port, "top_port")

        bottom_port = Port()
        bottom_port.transform(RotationAroundZ(pi))
        translate_to = map(operator.sub,
                first_part.bottom_port.middle.pos, bottom_port.middle.pos)
        bottom_port.transform(Translation(translate_to))
        self.add(bottom_port, "bottom_port")

if __name__ == "__main__":
    m = Pmpc(n=13)
    Plot(m).show()
