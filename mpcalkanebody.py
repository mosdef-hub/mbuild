__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *
from mpc import *

class MpcAlkaneBody(Compound):

    def __init__(self, ctx={}, direction='left'):
        super(MpcAlkaneBody, self).__init__(ctx=ctx)

        self.add(Atom(kind='CB', pos=(0, 0, 0)),'c')

        if direction == 'left':
            # mpc is on the left
            self.add(Atom(kind='HB', pos=(1, 0, 0)),'h1')
            self.add(Atom(kind='G', pos=(-1.2, 0, 0)),'c0') # this is where we want the MPC's first C atom to be
            self.add(Port(),'mpc_port')
            self.mpc_port.transform(RotationAroundY(pi/2))
            self.mpc_port.transform(Translation((-0.0,0,0)))
        else:
            # mpc is on the right
            self.add(Atom(kind='HB', pos=(-1, 0, 0)),'h1')
            self.add(Atom(kind='G', pos=(1.2, 0, 0)),'c0') # this is where we want the MPC's first C atom to be
            self.add(Port(),'mpc_port')
            self.mpc_port.transform(RotationAroundY(-pi/2))
            self.mpc_port.transform(Translation((0.0,0,0)))


        self.add(Mpc(), 'mpc')

        # transform using equivalence of two ports
        # m.mpc.transform([(m.mpc.mpc_port, m.mpc_port)])

        # transform using multiple equivalence relations (with a least squares solution)
        #      m.mpc.mpc_port should match up with m.mpc_port
        #    and
        #      the first C atom that was read from the xyz file should match up with our c0
        self.mpc.transform([(self.mpc.mpc_port, self.mpc_port),
                         (self.mpc.mpc_xyz.C_0, self.c0)])

        self.add(Port(),'male_port')
        self.male_port.transform(RotationAroundZ(pi))
        self.male_port.transform(Translation((0,-0.7,0)))

        self.add(Port(),'female_port')
        self.female_port.transform(RotationAroundZ(pi))
        self.female_port.transform(Translation((0,0.7,0)))


if __name__ == "__main__":
    m = MpcAlkaneBody.create()
    # m = Methane.create()
    # print ethane
    print m
    print m.atoms()
    # print m.label()
    m.plot(labels=False, verbose=False)
