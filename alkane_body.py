__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *
import pdb

class AlkaneBody(Compound):

    def __init__(self, ctx={}, direction='left'):
        super(AlkaneBody, self).__init__(ctx=ctx)
        self.add(Atom(kind='C', pos=(0, 0, 0)),'c')
        if direction == 'left':
            # self.add(Atom(kind='H', pos=(0.7, 0.0, 0.7)),'h1')
            # self.add(Atom(kind='H', pos=(-0.7, 0.0, 0.7)),'h2')
            self.add(Atom(kind='H', pos=(0.7, 0.0, 0.7)))
            self.add(Atom(kind='H', pos=(-0.7, 0.0, 0.7)))
        elif direction == 'right':
            # self.add(Atom(kind='H', pos=(0.7, 0, -0.7)),'h1')
            # self.add(Atom(kind='H', pos=(-0.7, 0, -0.7)),'h2')
            self.add(Atom(kind='H', pos=(0.7, 0, -0.7)))
            self.add(Atom(kind='H', pos=(-0.7, 0, -0.7)))
        else:
            raise Exception("Unknown 'direction' value: " + str(direction))

        self.add(Port(),'male_port')
        self.male_port.transform(RotationAroundZ(pi))
        self.male_port.transform(Translation((0,-0.7,0)))


        self.add(Port(),'female_port')
        self.female_port.transform(RotationAroundZ(pi))
        self.female_port.transform(Translation((0,0.7,0)))


