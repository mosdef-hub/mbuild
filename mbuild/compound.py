__author__ = 'sallai'
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import pdb
from mayavi import mlab
from tvtk.api import tvtk
import numpy as np
from mbuild.coordinate_transform import *
from mbuild.atom import *
from collections import OrderedDict
from itertools import *

class Compound(object):
    @classmethod
    def create(cls, label=None):
        """
        Create a compound
        :param label: a text label of the compount
        :return: the compound object
        """
        m = cls()
        m._components = OrderedDict() # this contain label to Compound or label to Atom mappings
        if label is not None:
            m._label = label
        return m

    def add(self, what, label=None):
        """
        Add a component to a Compound, inserting it into the map with a label as key, as well as adding it as a member
        variable with the label as variable name
        :param what: an Atom or Compound to be added
        :param label: optionally override the component's label
        :raise: Exception if component with the given label already exists in the map
        """
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
        Compute an equivalence transformation that transforms this compound to another compound's coordinate system.
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
                for label0, atom0 in pair[0].atoms():
                    atom1 = pair[1].component(label0)
                    self_points = vstack([self_points, atom0.pos])
                    other_points = vstack([other_points, atom1.pos])

        T = RigidTransform(self_points, other_points)
        return T

    def transform(self, T):
        """
        Transform this point cloud to another's coordinate system, or apply a given coordinate transformation.
        :param T: list of equivalence relations or coordinate transform
        """

        if not isinstance(T, CoordinateTransform):
            # we're assuming here that T is a list of equivalence relations
            T = self.createEquivalenceTransform(T)

        # # transform the contained atoms recursively
        # for label, component in self._components.iteritems():
        #     component.transform(T)

        # transform the contained atoms in batch
        arr = np.fromiter(chain.from_iterable(atom.pos for label, atom in self.atoms()), dtype=np.float64)
        arrnx3 = arr.reshape((-1,3))

        arrnx3 = T.apply(arrnx3)
        arr = arrnx3.reshape((-1))

        # write back new coordinates into atoms
        i = 0
        for label, atom in self.atoms():
            atom.pos = (arr[i], arr[i+1], arr[i+2])
            i=i+3

        return self

    # def atoms(self):
    #     """
    #     Get all atoms of the Compound recursively
    #     :return: map containing atoms with labels as keys
    #     """
    #     atoms = OrderedDict() # empty dict
    #     for label, component in self._components.iteritems():
    #         # add local atoms
    #         if isinstance(component, Atom):
    #             atoms[label] = component
    #             # add atoms in sub-components recursively
    #         if isinstance(component, Compound):
    #             for sublabel, subatom in component.atoms().iteritems():
    #                 atoms[label + '.' + sublabel] = subatom
    #     return atoms

    def atoms(self):
        """
        Generate all atoms of the Compound recursively
        :return: label - atom pairs
        """

        for label, component in self._components.iteritems():
            # add local atoms
            if isinstance(component, Atom):
                yield label, component
                # add atoms in sub-components recursively
            if isinstance(component, Compound):
                for sublabel, subatom in component.atoms():
                    yield label + '.' + sublabel, subatom

    def savexyz(self, fn, print_ports=False):
        """
        Save into an xyz file
        :param fn: file name
        :param print_ports: if False, ghost points are not written
        """
        with open(fn, 'w') as f:
            if print_ports:
                f.write(str(self.atoms().__len__()) + '\n\n')
            else:
                i = 0
                for key, value in self.atoms().iteritems():
                    if value.kind != 'G':
                        i += 1
                f.write(str(i) + '\n\n')
            for key, value in self.atoms().iteritems():
                if print_ports:
                    f.write(value.kind + '\t' +
                            str(value.pos[0]) + '\t' +
                            str(value.pos[1]) + '\t' +
                            str(value.pos[2]) + '\n')
                else:
                    if value.kind != 'G':
                        f.write(value.kind + '\t' +
                                str(value.pos[0]) + '\t' +
                                str(value.pos[1]) + '\t' +
                                str(value.pos[2]) + '\n')

    def component(self, component_path):
        """
        Find a component by label
        :param component_path: the label of the component (may be hierarchical)
        :return: the component (if found), None otherwise
        """
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
        """
        Compute the bounding box of the compound
        :rtype : (minx, miny, minz), (maxx, maxy, maxz)
        """
        minx = float('inf')
        miny = float('inf')
        minz = float('inf')
        maxx = float('-inf')
        maxy = float('-inf')
        maxz = float('-inf')

        for label, a in self.atoms():
            if excludeG and a.kind == 'G':
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
        """
        Plot atoms in 3d
        :param verbose: if True, atom types will be plotted
        :param labels: if True, labels will be plotted
        """

        x = []
        y = []
        z = []
        r = []
        g = []
        b = []
        rgb = []

        # sort atoms by type
        d = dict()

        for (label, atom) in self.atoms():
            if atom.kind != 'G' or verbose:
                if not atom.kind in d.keys():
                    d[atom.kind] = [atom]
                else:
                    d[atom.kind].append(atom);

        for (kind,atomList) in d.items():
            x = []
            y = []
            z = []
            r = []
            for atom in atomList:
                x.append(atom.pos[0])
                y.append(atom.pos[1])
                z.append(atom.pos[2])
                r.append(atom.vdw_radius)

            fig = mlab.points3d(x,y,z,r,color=atomList[0].colorRGB, scale_factor=1, scale_mode='scalar')
            #fig.glyph.glyph.clamping = False
        mlab.show()

    def plot2(self, verbose=False, labels=False):
        """
        Plot atoms in 3d
        :param verbose: if True, atom types will be plotted
        :param labels: if True, labels will be plotted
        """
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection='3d', aspect='equal')
        coord_min = inf
        coord_max = -inf
        for (label, atom) in self.atoms():
            if atom.kind != 'G' or verbose:
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
