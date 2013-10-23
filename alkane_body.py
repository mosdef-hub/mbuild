__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *

class AlkaneBody(Compound):

    @classmethod
    def create(cls, label=None):
        m = super(AlkaneBody, cls).create(label)

        m.add(C((0, 0, 0)),'c')
        m.add(H((1, 0, 0)),'h1')
        m.add(H((-1, 0, 0)),'h2')

        m.add(Port.create(),'male_port')
        m.male_port.applyTransformation(CoordinateTransform.rotation_around_z(pi))
        m.male_port.applyTransformation(CoordinateTransform.translation((0,-0.6,0)))


        m.add(Port.create(),'female_port')
        m.female_port.applyTransformation(CoordinateTransform.rotation_around_z(pi))
        m.female_port.applyTransformation(CoordinateTransform.translation((0,0.6,0)))
        return m


