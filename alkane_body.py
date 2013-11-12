__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *

class AlkaneBody(Compound):

    @classmethod
    def create(cls, ctx={}):
        m = super(AlkaneBody, cls).create(ctx=ctx)

        m.add(C((0, 0, 0)),'c')
        m.add(H((1, 0, 0)),'h1')
        m.add(H((-1, 0, 0)),'h2')

        m.add(Port.create(),'male_port')
        m.male_port.transform(RotationAroundZ(pi))
        m.male_port.transform(Translation((0,-0.7,0)))


        m.add(Port.create(),'female_port')
        m.female_port.transform(RotationAroundZ(pi))
        m.female_port.transform(Translation((0,0.7,0)))
        return m


