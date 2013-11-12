__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *
from mpc import *

class MpcAlkaneBody(Compound):

    @classmethod
    def create(cls, label=None, direction='left', ctx={}):
        m = super(MpcAlkaneBody, cls).create(label)

        m.add(C((0, 0, 0)),'c')

        if direction == 'left':
            # mpc is on the left
            m.add(H((1, 0, 0)),'h1')
            m.add(G((-1.2, 0, 0)),'c0') # this is where we want the MPC's first C atom to be
            m.add(Port.create(),'mpc_port')
            m.mpc_port.transform(RotationAroundY(pi/2))
            m.mpc_port.transform(Translation((-0.0,0,0)))
        else:
            # mpc is on the right
            m.add(H((-1, 0, 0)),'h1')
            m.add(G((1.2, 0, 0)),'c0') # this is where we want the MPC's first C atom to be
            m.add(Port.create(),'mpc_port')
            m.mpc_port.transform(RotationAroundY(-pi/2))
            m.mpc_port.transform(Translation((0.0,0,0)))


        m.add(Mpc.create(), 'mpc')

        # transform using equivalence of two ports
        # m.mpc.transform([(m.mpc.mpc_port, m.mpc_port)])

        # transform using multiple equivalence relations (with a least squares solution)
        #      m.mpc.mpc_port should match up with m.mpc_port
        #    and
        #      the first C atom that was read from the xyz file should match up with our c0
        m.mpc.transform([(m.mpc.mpc_port, m.mpc_port),
                         (m.mpc.mpc_xyz.C0, m.c0)])

        m.add(Port.create(),'male_port')
        m.male_port.transform(RotationAroundZ(pi))
        m.male_port.transform(Translation((0,-0.7,0)))

        m.add(Port.create(),'female_port')
        m.female_port.transform(RotationAroundZ(pi))
        m.female_port.transform(Translation((0,0.7,0)))
        return m


if __name__ == "__main__":
    m = MpcAlkaneBody.create()
    # m = Methane.create()
    # print ethane
    print m
    print m.atoms()
    # print m.label()
    m.plot(labels=False, verbose=False)
