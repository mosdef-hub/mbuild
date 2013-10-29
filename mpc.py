__author__ = 'sallai'
from mbuild.compound import *
from mbuild.xyz import *
from mbuild.port import *

class Mpc(Compound):

    @classmethod
    def create(cls, label=None):
        m = super(Mpc, cls).create(label)
        mpc = Xyz.create('mpc.xyz')
        m.add(mpc,'mpc_xyz')

        m.add(Port.create(),'mpc_port')
        #m.mpc_port.applyTransformation(CoordinateTransform.rotation_around_z(pi))
        #m.mpc_port.applyTransformation(CoordinateTransform.rotation_around_y(pi/2))
        m.mpc_port.applyTransformation(CoordinateTransform.translation((-2.2,-2.49,1.17)))

        return m

if __name__ == "__main__":
    m = Mpc.create(label='mpc')
    # m = Methane.create()
    # print ethane
    print m
    print m.atoms()
    print m.label()
    m.plot(labels=False, verbose=False)
