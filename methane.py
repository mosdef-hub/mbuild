__author__ = 'sallai'
from mbuild.compound import *


class Methane(Compound):

    @classmethod
    def create(cls, label=None):
        m = super(Methane, cls).create(label)
        m.add(C((0, 0, 0)),'c')
        m.add(H((1, 0, 0)),'h1')
        m.add(H((0, 1, 0)),'h2')
        m.add(H((-1, 0, 0)),'h3')
        m.add(H((0, -1, 0)),'h4')
        return m

if __name__ == "__main__":
    m = Methane.create(label='myMethane')
    # m = Methane.create()
    # print ethane
    print m
    print m.atoms()
    print m.label()
    m.plot()
