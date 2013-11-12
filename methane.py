__author__ = 'sallai'
from mbuild.compound import *


class Methane(Compound):

    @classmethod
    def create(cls, ctx={}):
        m = super(Methane, cls).create(ctx=ctx)
        m.add(C((0, 0, 0)),'c')
        m.add(H((1, 0, 0)),'h1')
        m.add(H((0, 1, 0)),'h2')
        m.add(H((-1, 0, 0)),'h3')
        m.add(H((0, -1, 0)),'h4')
        return m

if __name__ == "__main__":
    self = Methane.create()
    # m = Methane.create()
    # print ethane
    print self
    print self.atoms()
    # print m.label()
    self.plot()
