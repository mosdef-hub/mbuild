from math import sqrt

from mbuild.coordinate_transform import *
from mbuild.atom import Atom
from mbuild.compound import Compound
from mbuild.port import Port

class AlkaneTail(Compound):

    def __init__(self, ctx={}, ff='opls'):
        super(AlkaneTail, self).__init__(ctx=ctx)
        s2 = sqrt(2) / 2.0
        if ff == 'opls':
            h_kind = 'opls_140'
            c_kind = 'opls_135'
            # TODO: need differentiation depending on bonding structure
        else:
            h_kind = 'H'
            c_kind = 'C'

        self.add(Atom(kind=h_kind, pos=(s2, s2, -s2)), 'h1')
        self.add(Atom(kind=h_kind, pos=(s2, s2, s2)), 'h2')
        self.add(Atom(kind=h_kind, pos=(-s2, s2, s2)), 'h3')
        self.add(Atom(kind=c_kind, pos=(0, 0, 0)), 'c')

        self.add(Port(), 'female_port')
        self.female_port.transform(Translation((0,-0.7,0)))

if __name__ == "__main__":
    m = AlkaneTail()
    from mbuild.plot import Plot
    Plot(m, verbose=True).show()
