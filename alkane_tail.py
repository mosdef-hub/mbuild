__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *
from mbuild.plot import Plot

class AlkaneTail(Compound):

    def __init__(self, ctx={}):
        super(AlkaneTail, self).__init__(ctx=ctx)
        self.add(Atom(kind='H', pos=(1, 0, 0)), 'h1')
        self.add(Atom(kind='H', pos=(0, 1, 0)), 'h2')
        self.add(Atom(kind='H', pos=(-1, 0, 0)), 'h3')
        self.add(Atom(kind='C', pos=(0, 0, 0)), 'c')

        self.add(Port(), 'bottom_port')
        self.bottom_port.transform(Translation((0,-0.7,0)))


if __name__ == "__main__":
    m = AlkaneTail()
    # TreeView(m, verbose=True).show()
    Plot(m, verbose=True).show()
