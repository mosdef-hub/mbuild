__author__ = 'sallai'

import inspect
import sys

from mbuild.coordinate_transform import *

class Atom(object):

    def __init__(self, kind, pos=(0, 0, 0), charge=0):
        assert(isinstance(kind, basestring))

        self.kind = kind
        self.pos = pos
        self.charge = charge
        self.bonds = set()
        self.angles = set()
        self.dihedrals = set()

    @staticmethod
    def distance(a1, a2):
        return a1.distance(a2)

    def distance(self, a2):
        return sqrt((self.pos[0]-a2.pos[0])**2 + (self.pos[1]-a2.pos[1])**2 + (self.pos[2]-a2.pos[2])**2)

    def transform(self, T):
        self.pos = tuple(squeeze(T.apply(array([self.pos]))))

    def __repr__(self):
        return "Atom"+str(id(self))+"(" + self.kind + "," + str(self.pos) + ")"