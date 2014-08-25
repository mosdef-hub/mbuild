__author__ = 'sallai'
from mbuild.file_formats.xyz import *
from mbuild.port import *

class Mpc(Compound):

    def __init__(self, ctx={}):
        super(Mpc, self).__init__(ctx=ctx)

        mpc = Xyz('mpc.xyz')
        self.add(mpc,'mpc_xyz')

        self.add(Port(),'mpc_port')
        self.mpc_port.transform(Translation((-2.2,-2.49,1.17)))

if __name__ == "__main__":
    m = Mpc()
    # m = Methane.create()
    # print ethane
    print m
    # print m.label()
    TreeView(m).show()
