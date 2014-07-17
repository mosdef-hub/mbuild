import operator
import pdb

from numpy import pi
from mbuild.bond import Bond

from mbuild.coordinate_transform import *
from mbuild.compound import Compound
from mbuild.port import Port

from mpc_monomer import MpcMonomer


class Pmpc(Compound):

    def __init__(self, n=2, alpha=pi/6):
        if n < 1:
            raise Exception('n must be 1 or more')
        Compound.__init__(self)

        # n-1 times the body CH_2
        last_part = None
        for body_count in range(0, n):
            this_part = MpcMonomer(alpha=alpha)
            self.add(this_part, 'monomer_{0}'.format(body_count))
            if last_part is None:
                first_part = this_part
            else:
                transform(this_part, [(this_part.labels['bottom_port'],
                                      last_part.labels['top_port'])])
                bond = Bond(this_part.labels['C.3_2'], last_part.labels['C.3_38'])
                self.add(bond)


            last_part = this_part

        top_port = Port()
        transform(top_port, RotationAroundZ(pi))

        translate_to = last_part.labels['top_port'].labels['middle'].pos - top_port.labels['middle'].pos


        transform(top_port, Translation(translate_to))
        self.add(top_port, "top_port")

        bottom_port = Port()
        transform(bottom_port, RotationAroundZ(pi))

        translate_to = first_part.labels['bottom_port'].labels['middle'].pos - bottom_port.labels['middle'].pos
        transform(bottom_port, Translation(translate_to))
        self.add(bottom_port, "bottom_port")



if __name__ == "__main__":
    m = Pmpc(n=13)
    from mbuild.plot import Plot
    Plot(m).show()
