__author__ = 'sallai'

import inspect
import sys

from mbuild.coordinate_transform import *

class Atom(object):

    # subclasses = dict()
    vdw_radius=1.0
    color='w'
    colorRGB=(1,1,1)

    def __init__(self, kind, pos=(0, 0, 0)):

        if not isinstance(kind, basestring):
            raise Exception("Kind must be as a non-empty string")

        self.pos = pos
        self.bonds = set()
        self.angles = set()
        self.dihedrals = set()

        self.kind = kind

        # # create a class for this kind, add it to the subclasses dict (kind to class mappings)
        # if not self.__class__.subclasses.has_key(kind):
        #     self.__class__.subclasses[kind] = type(kind, (Atom,), {"kind": kind})
        #
        # self.__class__ = self.__class__.subclasses[kind]

    # def setKind(self, kind):
    #     # if it doesn't exist, create a class for this kind, add it to the subclasses dict (kind to class mappings)
    #     if not self.__class__.subclasses.has_key(kind):
    #         self.__class__.subclasses[kind] = type(kind, (Atom,), {"kind": kind})
    #
    #     # set the class of the atom to the class of the kind
    #     self.__class__ = self.__class__.subclasses[kind]


    # def defaultLabel(self):
    #     return self.kind + '_' + str(id(self))
    #
    # @classmethod
    # def create(cls, kind=None, pos=(0, 0, 0)):
    #     # try instantiating a dynamically created class
    #     for classname, classobj in inspect.getmembers(sys.modules[cls.__module__], inspect.isclass):
    #         if classname == kind:
    #             return classobj(pos=pos)
    #     # resort to Atom if no specialized class exists for the given kind
    #     return Atom(kind=kind, pos=pos)

    @staticmethod
    def distance(a1, a2):
        return a1.distance(a2)

    def distance(self, a2):
        return sqrt((self.pos[0]-a2.pos[0])**2 + (self.pos[1]-a2.pos[1])**2 + (self.pos[2]-a2.pos[2])**2)

    def transform(self, T):
        self.pos = tuple(squeeze(T.apply(array([self.pos]))))

    def __repr__(self):
        return "Atom"+str(id(self))+"(" + self.kind + "," + str(self.pos) + ")"


# def AtomClassFactory(name, vdw_radius=1.0, color='w', colorRGB=(1,1,1)):
#     # define a prototype init function for the new class that calls the init function of Atom
#     def __init__(self, pos=None):
#         Atom.__init__(self, kind=name, pos=pos)
#     # create new class, setting its init function to the prototype
#     newclass = type(name, (Atom,),{"__init__": __init__, "kind": name, "vdw_radius": vdw_radius, "color": color, "colorRGB": colorRGB})
#     return newclass
#
# class G(Atom):
#     vdw_radius = 1.2
#     color = "gray"
#     colorRGB = (.5, .5, .5)
#
#     def __init__(self, pos=(0, 0, 0)):
#         super(G, self).__init__('G', pos)
#
#
# C = AtomClassFactory('C', vdw_radius=1.7, color='teal', colorRGB=(0, .5, .5))
# CB = AtomClassFactory('CB', vdw_radius=1.7, color='teal', colorRGB=(0, .5, .5))
# C2 = AtomClassFactory('C2', vdw_radius=1.7, color='teal', colorRGB=(0, 0, 0))
# H = AtomClassFactory('H', vdw_radius=1.2, color='white', colorRGB=(1, 1, 1))
# HB = AtomClassFactory('HB', vdw_radius=1.2, color='white', colorRGB=(1, 1, 1))
# O = AtomClassFactory('O', vdw_radius=1.52, color='red', colorRGB=(1, 0, 0))
# O1 = AtomClassFactory('O', vdw_radius=1.52, color='red', colorRGB=(1, 0, 0))
# F = AtomClassFactory('F', vdw_radius=1.35, color='pink', colorRGB=(1, 192.0/256.0, 203.0/256.0))
# N = AtomClassFactory('N', vdw_radius=1.55, color='blue', colorRGB=(0, 0, 1))
# P = AtomClassFactory('P', vdw_radius=1.8, color='orange', colorRGB=(1, 165.0/256.0, 0))
# Si = AtomClassFactory('Si', vdw_radius=2.10, color='yellow', colorRGB=(1, 1, 0))
# Si1 = AtomClassFactory('Si', vdw_radius=2.10, color='yellow', colorRGB=(1, 1, 0))
