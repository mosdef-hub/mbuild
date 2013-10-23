__author__ = 'sallai'
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

from mbuild.coordinate_transform import *

from mbuild.atom import *


class Compound(object):
    @classmethod
    def create(cls, label=None):
        m = cls()
        m._components = {}
        if label is not None:
            m._label = label
        return m

    def add(self, what, label=None):
        if label is None:
            label = what.label()
        if label in self._components.keys():
            raise Exception("label " + label + " already exists in " + str(what))
        self._components[label] = what
        setattr(self, label, what)

    def label(self):
        if hasattr(self, '_label'):
            return self._label
        else:
            return str(self.__class__.__name__) + '_' + str(id(self))

    def createEquivalenceTransform(self, equiv):
        """
        Compute an equivalence transformation that transforms this point cloud to another's coordinate system.
        :param other: the other point cloud
        :param equiv: list of equivalent points
        :returns: the coordinatetransform object that transforms this point cloud to the other point cloud's coordinates system
        """

        self_points = array([])
        self_points.shape = (0, 3)
        other_points = array([])
        other_points.shape = (0, 3)

        for pair in equiv:
            if not isinstance(pair, tuple) or len(pair) != 2:
                raise Exception('Equivalence pair not a 2-tuple')
            if not ((isinstance(pair[0], Compound) and isinstance(pair[1], Compound)) or (
                    isinstance(pair[0], Atom) and isinstance(pair[1], Atom))):
                raise Exception(
                    'Equivalence pair type mismatch: pair[0] is a ' + str(type(pair[0])) + ' and pair[1] is a ' + str(
                        type(pair[1])))

            if isinstance(pair[0], Atom):
                self_points = vstack([self_points, pair[0].pos])
                other_points = vstack([other_points, pair[1].pos])
            if isinstance(pair[0], Compound):
                for label0, atom0 in pair[0].atoms().iteritems():
                    atom1 = pair[1].component(label0)
                    self_points = vstack([self_points, atom0.pos])
                    other_points = vstack([other_points, atom1.pos])

        T = CoordinateTransform.compute(self_points, other_points)
        return T

    def transform(self, equiv):
        """
        Transform this point cloud to another's coordinate system.
        :param other: the other point cloud
        :param equiv: list of equivalent points
        :returns: the matrix that transforms this point cloud to the other point cloud's coordinates system
        """
        T = self.createEquivalenceTransform(equiv)

        # transform the contained atoms recursively
        self.applyTransformation(T)

        return self

    def applyTransformation(self, T):
        for label, component in self._components.iteritems():
            component.applyTransformation(T)


    def atoms(self):
        atoms = {} # empty dict
        for label, component in self._components.iteritems():
            # add local atoms
            if isinstance(component, Atom):
                atoms[label] = component
                # add atoms in sub-components recursively
            if isinstance(component, Compound):
                for sublabel, subatom in component.atoms().iteritems():
                    atoms[label + '.' + sublabel] = subatom
        return atoms

    def savexyz(self, fn):
        with open(fn, 'w') as f:
            f.write(str(self.atoms().__len__()) + '\n\n')
            for key, value in self.atoms().iteritems():
                f.write(value.atomType + '\t' + str(value.pos[0]) + '\t' + str(value.pos[1]) + '\t' + str(
                    value.pos[2]) + '\n')


    def component(self, component_path):
        dot_pos = component_path.find('.')
        if dot_pos > -1:
            subcomponent_path = component_path[:dot_pos]
            subpath = component_path[dot_pos + 1:]
            if subcomponent_path in self._components.keys():
                return self._components[subcomponent_path].component(subpath)
            else:
                return None
        else:
            if component_path in self._components.keys():
                return self._components[component_path]
            else:
                return None

    def boundingbox(self, excludeG=True):
        minx = float('inf')
        miny = float('inf')
        minz = float('inf')
        maxx = float('-inf')
        maxy = float('-inf')
        maxz = float('-inf')

        for label, a in self.atoms().iteritems():
            if excludeG and a.atomType == 'G':
                continue
            if a.pos[0] < minx:
                minx = a.pos[0]
            if a.pos[0] > maxx:
                maxx = a.pos[0]
            if a.pos[1] < miny:
                miny = a.pos[1]
            if a.pos[1] > maxy:
                maxy = a.pos[1]
            if a.pos[2] < minz:
                minz = a.pos[2]
            if a.pos[2] > maxz:
                maxz = a.pos[2]

        return (minx, miny, minz), (maxx, maxy, maxz)

    def plot(self, verbose=False, labels=True):
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection='3d', aspect='equal')
        for (label, atom) in self.atoms().items():
            if atom.atomType != 'G' or verbose:
                # print atom
                if labels:
                    atom.plot(ax, str(atom))
                else:
                    atom.plot(ax, None)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(self.label())

        pyplot.show()