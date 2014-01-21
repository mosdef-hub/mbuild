__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *
import pdb

class AlkaneBody(Compound):

    @classmethod
    def create(cls, ctx={}, direction='left'):
        m = super(AlkaneBody, cls).create(ctx=ctx)

        m.add(CB((0, 0, 0)),'c')
        if direction == 'left':
            m.add(HB((0.7, 0.0, 0.7)),'h1')
            m.add(HB((-0.7, 0.0, 0.7)),'h2')
        elif direction == 'right':
            m.add(HB((0.7, 0, -0.7)),'h1')
            m.add(HB((-0.7, 0, -0.7)),'h2')
        else:
            raise Exception("Unknown 'direction' value: " + str(direction))

        m.add(Port.create(),'male_port')
        m.male_port.transform(RotationAroundZ(pi))
        m.male_port.transform(Translation((0,-0.7,0)))


        m.add(Port.create(),'female_port')
        m.female_port.transform(RotationAroundZ(pi))
        m.female_port.transform(Translation((0,0.7,0)))
        return m


