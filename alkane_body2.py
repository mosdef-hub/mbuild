__author__ = 'sallai'
from math import pi
from mbuild.coordinate_transform import RotationAroundZ, Translation
from mbuild.compound import Compound
from mbuild.port import Port
import mbuild.atom as a

class alkane_body2(Compound):

    @classmethod
    def create(cls, ctx={}, cwd="", n=1):
        m = super(alkane_body2, cls).create(ctx=ctx)

        print n

        m.add(a.C((0, 0, 0)),'c')
        m.add(a.H((1, 0, 0)),'h1')
        m.add(a.H((-1, 0, 0)),'h2')

        m.add(Port.create(),'male_port')
        m.male_port.transform(RotationAroundZ(pi))
        m.male_port.transform(Translation((0,-0.7,0)))


        m.add(Port.create(),'female_port')
        m.female_port.transform(RotationAroundZ(pi))
        m.female_port.transform(Translation((0,0.7,0)))
        return m


print "alkane_body2.py loaded"