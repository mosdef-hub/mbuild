__author__ = 'sallai'
from mbuild.coordinate_transform import *


class Atom(object):
    def __init__(self, atomType=None, pos=(0, 0, 0)):
        super(Atom, self).__init__()
        self.atomType = atomType
        self.pos = pos
        self.color = 'k'

    def label(self):
        return self.atomType + '_' + str(id(self))

    @classmethod
    def create(cls, atomType=None, pos=(0, 0, 0)):
        if atomType == 'C':
            return C(pos)
        elif atomType == 'H':
            return H(pos)
        elif atomType == 'O':
            return O(pos)
        elif atomType == 'F':
            return F(pos)
        else:
            return Atom(atomType, pos)


    def distance(self, a2):
        return sqrt((self.pos[0]-a2.pos[0])**2 + (self.pos[1]-a2.pos[1])**2 + (self.pos[2]-a2.pos[2])**2)

    def applyTransformation(self, T):
        self.pos = tuple(squeeze(T.transform(array([self.pos]))))

    def __repr__(self):
        return "Atom[(" + self.atomType + ")" + str(self.pos) + "]"

    def plot(self, ax, text):
        # print atom
        ax.scatter(self.pos[0], self.pos[1], self.pos[2], c=self.color, marker='o')
        if text is not None:
            ax.text(self.pos[0], self.pos[1], self.pos[2],  text)
        else:
            ax.text(self.pos[0], self.pos[1], self.pos[2],  self.atomType)


class C(Atom):
    def __init__(self, pos=(0, 0, 0)):
        super(C, self).__init__('C', pos)
        self.color = 'y'


class H(Atom):
    def __init__(self, pos=(0, 0, 0)):
        super(H, self).__init__('H', pos)
        self.color = 'w'


class O(Atom):
    def __init__(self, pos=(0, 0, 0)):
        super(O, self).__init__('O', pos)
        self.color = 'b'

class F(Atom):
    def __init__(self, pos=(0, 0, 0)):
        super(F, self).__init__('F', pos)
        self.color = 'g'


class N(Atom):
    def __init__(self, pos=(0, 0, 0)):
        super(N, self).__init__('N', pos)


class G(Atom):
    def __init__(self, pos=(0, 0, 0)):
        super(G, self).__init__('G', pos)

