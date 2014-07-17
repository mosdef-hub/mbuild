from copy import deepcopy
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

        proto = MpcMonomer(alpha=alpha)

        last_part = None
        for body_count in range(0, n):
            this_part = deepcopy(proto)
            self.add(this_part, 'monomer_{0}'.format(body_count))
            if last_part is None:
                first_part = this_part
            else:
                # transform this part, such that it's bottom port is rotated+translated to the last part's top port
                equivalence_transform(this_part, this_part.bottom_port,
                                      last_part.top_port)
                bond = Bond(this_part.labels['C.3_2'], last_part.labels['C.3_38'])
                self.add(bond)
            last_part = this_part

        # hoist the last part's top port to be the top port of the polymer
        self.add(last_part.top_port, 'top_port', containment=False)

        # hoist the first part's bottom port to be the bottom port of the polymer
        self.add(first_part.bottom_port, 'bottom_port', containment=False)

if __name__ == "__main__":
    m = Pmpc(n=13)
    from mbuild.plot import Plot
    Plot(m).show()
