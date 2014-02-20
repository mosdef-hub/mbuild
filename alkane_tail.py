__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *

class AlkaneTail(Compound):

    def __init__(self, ctx={}):
        super(AlkaneTail, self).__init__(ctx=ctx)
        self.add(Atom(kind='H', pos=(1, 0, 0)), 'h1')
        self.add(Atom(kind='H', pos=(0, 1, 0)), 'h2')
        self.add(Atom(kind='H', pos=(-1, 0, 0)), 'h3')
        self.add(Atom(kind='C', pos=(0, 0, 0)), 'c')

        self.add(Port(), 'female_port')
        self.female_port.transform(Translation((0,-0.7,0)))

        self.add(Port(), 'male_port')
        self.male_port.transform(RotationAroundZ(pi))
        self.male_port.transform(Translation((0,-0.7,0)))

