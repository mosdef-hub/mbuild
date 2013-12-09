__author__ = 'sallai'

import inspect
import sys

from mbuild.coordinate_transform import *

class Atom(object):

    vdw_radius=1.0
    color='w'
    colorRGB=(1,1,1)

    def __init__(self, kind=None, pos=(0, 0, 0)):
        super(Atom, self).__init__()
        self.kind = kind
        self.pos = pos

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

    # def __hash__(self):
    #     h = hash((self.kind, self.pos))
    #     return (h ^ (h << 16) ^ 89869747) * 3644798167
    #
    # def __eq__(self, other):
    #     return self.__hash__() == other.__hash__()

    def __repr__(self):
        return "Atom[(" + self.kind + ")" + str(self.pos) + "]"

    # def plot(self, ax, text):
    #     # print atom
    #     ax.scatter(self.pos[0], self.pos[1], self.pos[2],
    #             c=self.color,
    #             marker='o',
    #             s=100 * self.vdw_radius ** 3)
    #     if text is not None:
    #         ax.text(self.pos[0], self.pos[1], self.pos[2],  text)
    #     else:
    #         ax.text(self.pos[0], self.pos[1], self.pos[2],  self.kind)

    # def plot3_worker(self, g, vnames, parent_vi):
    #     vi = g.add_vertex()
    #     if parent_vi >= 0:
    #         g.add_edges((parent_vi,vi))

def AtomClassFactory(name, vdw_radius=1.0, color='w', colorRGB=(1,1,1)):
    # define a prototype init function for the new class that calls the init function of Atom
    def __init__(self, pos=None):
        Atom.__init__(self, kind=name, pos=pos)
    # create new class, setting its init function to the prototype
    newclass = type(name, (Atom,),{"__init__": __init__, "kind": name, "vdw_radius": vdw_radius, "color": color, "colorRGB": colorRGB})
    return newclass

class G(Atom):
    vdw_radius = .2
    color = "gray"
    colorRGB = (.5, .5, .5)

    def __init__(self, pos=(0, 0, 0)):
        super(G, self).__init__('G', pos)


C = AtomClassFactory('C', vdw_radius=1.7, color='teal', colorRGB=(0, .5, .5))
CB = AtomClassFactory('CB', vdw_radius=1.7, color='teal', colorRGB=(0, .5, .5))
C2 = AtomClassFactory('C2', vdw_radius=1.7, color='teal', colorRGB=(0, 0, 0))
H = AtomClassFactory('H', vdw_radius=1.2, color='white', colorRGB=(1, 1, 1))
HB = AtomClassFactory('HB', vdw_radius=1.2, color='white', colorRGB=(1, 1, 1))
O = AtomClassFactory('O', vdw_radius=1.52, color='red', colorRGB=(1, 0, 0))
F = AtomClassFactory('F', vdw_radius=1.35, color='pink', colorRGB=(1, 192.0/256.0, 203.0/256.0))
N = AtomClassFactory('N', vdw_radius=1.55, color='blue', colorRGB=(0, 0, 1))
P = AtomClassFactory('P', vdw_radius=1.8, color='orange', colorRGB=(1, 165.0/256.0, 0))
Si = AtomClassFactory('Si', vdw_radius=2.10, color='purple', colorRGB=(1, 0, 165.0/256.0))
