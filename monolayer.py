from mbuild.port import Port

__author__ = 'sallai'

from alkane_tail import AlkaneTail
from alkane_body import AlkaneBody
from mbuild.compound import *
from methane import *
from surface import Surface
from n_alkyl import NAlkyl

class Monolayer(Compound):
    @classmethod
    def create(cls, chain_length=4, coverage=1, ctx={}):

        m = super(Monolayer, cls).create(ctx=ctx)


        m.add(Surface.create(ctx=ctx),'surface')

        for label, port in m.surface._components.iteritems():
            if isinstance(port, Port):
                alkyl = NAlkyl.create(ctx=ctx, n=chain_length)
                alkyl.transform([(alkyl.port, port)])
                m.add( alkyl,'alkyl_for_' + label)

        return m

if __name__ == "__main__":
    m = Monolayer.create(chain_length=8)
    # print [(label,atom.pos) for label, atom in m.atoms()]
    m.plot(labels=False, verbose=True)