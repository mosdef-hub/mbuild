from mbuild.atom import Atom

__author__ = 'sallai'
from math import pi
from mbuild.coordinate_transform import RotationAroundZ, Translation
from mbuild.compound import Compound
from mbuild.port import Port
import mbuild.atom as a

class AlkaneBody(Compound):

    def __init__(self, ctx={}, cwd="", n=1):
        super(AlkaneBody, self).__init__(ctx=ctx)

        self.add(Atom(kind='C', pos=(0, 0, 0)),'c')
        self.add(Atom(kind='H', pos=(1.0, 0.0, 0.0)), 'h1')
        self.add(Atom(kind='H', pos=(-1.0, 0.0, 0.0)), 'h2')

        self.add(Port(),'male_port')
        self.male_port.transform(RotationAroundZ(pi))
        self.male_port.transform(Translation((0,-0.7,0)))


        self.add(Port(),'female_port')
        self.female_port.transform(RotationAroundZ(pi))
        self.female_port.transform(Translation((0,0.7,0)))
