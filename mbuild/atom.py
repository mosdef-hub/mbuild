__author__ = 'sallai'
from mbuild.coordinate_transform import *
import inspect
import sys

class Atom(object):
    def __init__(self, kind=None, pos=(0, 0, 0), wdv_radius=1.0, color='w', colorRGB=(1,1,1)):
        super(Atom, self).__init__()
        self.kind = kind
        self.pos = pos
        self.color = color
        self.colorRGB = colorRGB
        self.vdw_radius = wdv_radius

    def defaultLabel(self):
        return self.kind + '_' + str(id(self))

    @classmethod
    def create(cls, kind=None, pos=(0, 0, 0)):
        # try instantiating a dynamically created class
        for classname, classobj in inspect.getmembers(sys.modules[cls.__module__], inspect.isclass):
            if classname == kind:
                return classobj(pos=pos)
        # resort to Atom if no specialized class exists for the given kind
        return Atom(kind=kind, pos=pos)

    @staticmethod
    def distance(a1, a2):
        return a1.distance(a2)

    def distance(self, a2):
        return sqrt((self.pos[0]-a2.pos[0])**2 + (self.pos[1]-a2.pos[1])**2 + (self.pos[2]-a2.pos[2])**2)

    def transform(self, T):
        self.pos = tuple(squeeze(T.apply(array([self.pos]))))

    def __repr__(self):
        return "Atom[(" + self.kind + ")" + str(self.pos) + "]"

    def plot(self, ax, text):
        # print atom
        ax.scatter(self.pos[0], self.pos[1], self.pos[2],
                c=self.color,
                marker='o',
                s=100 * self.vdw_radius ** 3)
        if text is not None:
            ax.text(self.pos[0], self.pos[1], self.pos[2],  text)
        else:
            ax.text(self.pos[0], self.pos[1], self.pos[2],  self.kind)


def AtomClassFactory(name, wdv_radius=1.0, color='w', colorRGB=(1,1,1)):
    # define a prototype init function for the new class that calls the init function of Atom
    def __init__(self, pos=None):
        Atom.__init__(self, kind=name, pos=pos, wdv_radius=wdv_radius, color=color, colorRGB=colorRGB)
    # create new class, setting its init function to the prototype
    newclass = type(name, (Atom,),{"__init__": __init__})
    return newclass

class G(Atom):
    def __init__(self, pos=(0, 0, 0)):
        super(G, self).__init__('G', pos)
        self.colorRGB = (.5, .5, .5)

C = AtomClassFactory('C', wdv_radius=1.7, color='teal', colorRGB=(0, .5, .5))
H = AtomClassFactory('H', wdv_radius=1.2, color='white', colorRGB=(1, 1, 1))
O = AtomClassFactory('O', wdv_radius=1.52, color='red', colorRGB=(1, 0, 0))
F = AtomClassFactory('F', wdv_radius=1.35, color='pink', colorRGB=(1, 192/256, 203/256))
N = AtomClassFactory('N', wdv_radius=1.55, color='blue', colorRGB=(0, 0, 1))
P = AtomClassFactory('P', wdv_radius=1.8, color='orange', colorRGB=(1, 165/256, 0))

