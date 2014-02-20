from mbuild.port import Port
from mbuild.treeview import TreeView

__author__ = 'sallai'

from mbuild.compound import *
from methane import *
from surface import Surface
from n_alkyl import NAlkyl

class Monolayer(Compound):

    def __init__(self, chain_length=4, coverage=1, ctx={}):

        super(Monolayer, self).__init__(ctx=ctx)

        self.add(Surface(ctx=ctx),'surface')

        for port in self.surface.parts:
            if isinstance(port, Port):
                alkyl = NAlkyl(ctx=ctx, n=chain_length)
                alkyl.transform([(alkyl.port, port)])
                self.add( alkyl)

if __name__ == "__main__":
    m = Monolayer(chain_length=8)
    # print [(label,atom.pos) for label, atom in m.atoms()]
    TreeView(m).show()
    # m.plot(labels=False, verbose=True)