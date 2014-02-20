
from mbuild.angle import Angle
from mbuild.bond import Bond
from mbuild.dihedral import Dihedral
from mbuild.periodic_kdtree import PeriodicCKDTree

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
from copy import copy, deepcopy
from mbuild.orderedset import *
import numpy as np
from numpy.linalg import norm

class Compound(object):

    def __init__(self, ctx={}, kind=None, bounds = [0.0, 0.0, 0.0]):
        # set the context (optional)
        self.ctx = ctx
        self.bounds = bounds
        # set the kind (defaults to classname)
        if kind:
            self.kind = kind
        else:
            self.kind = self.__class__.__name__

        self.parts = OrderedSet() # contains children (compounds or atoms)
        self.references = OrderedDict()  # label to compound/atom mappings -- need not be in parts

        self.bonds = OrderedSet()
        self.angles = OrderedSet()
        self.dihedrals = OrderedSet()

    # def findLabelOfComponent(self, what):
    #     return self._components.keys()[self._components.values().index(what)]

    # def addAlias(self, what, alias_label):
    #     if isinstance(what, basestring):
    #         label = what
    #         if not label:
    #             raise Exception("no label specified")
    #     else:
    #         label = self.findLabelOfComponent(what)
    #         if not label:
    #             raise Exception("no label found for component")
    #
    #     if label != alias_label and alias_label in self._components.keys():
    #         raise Exception("label " + label + " already exists in " + str(self))
    #
    #
    #     self._aliases[alias_label] = label
    #
    #     setattr(self, alias_label, self._components[label])

    def add(self, new_obj, label=None):

        # add atom as a part, set a reference to it
        if isinstance(new_obj, Atom):
            self.parts.add(new_obj)
            if label is not None:

                # support label with counter
                if label.endswith("#"):
                    i = 0
                    while label.replace("#", str(i)) in self.references.keys():
                        i += 1
                    label = label.replace("#", str(i))

                if label in self.references.keys():
                    raise Exception("label " + label + " already exists in " + str(self))
                else:
                    self.references[label] = new_obj
                    setattr(self, label, new_obj)

        # add compound as a part, set a reference to it
        elif isinstance(new_obj, Compound):
            self.parts.add(new_obj)
            if label is not None:

                # support label with counter
                if label.endswith("#"):
                    i = 0
                    while label.replace("#", str(i)) in self.references.keys():
                        i += 1
                    label = label.replace("#", str(i))

                if label in self.references.keys():
                    raise Exception("label " + label + " already exists in " + str(self))
                else:
                    self.references[label] = new_obj
                    setattr(self, label, new_obj)

        # add bond to compound
        elif isinstance(new_obj, Bond):
            # don't add B-A if A-B is already added
            # print "self.bonds before adding: " + str(list(self.bonds))
            # print "adding bond "+str(new_obj)
            if not new_obj in self.bonds and not new_obj.cloneImage() in self.bonds:
                self.bonds.add(new_obj)
                new_obj.atom1.bonds.add(new_obj)
                new_obj.atom2.bonds.add(new_obj)
            # print "self.bonds after adding:  " + str(list(self.bonds))


        # add angle to compound
        elif isinstance(new_obj, Angle):
            # don't add A-B-C if C-B-A is already added
            if not new_obj in self.angles and not new_obj.cloneImage() in self.angles:
                self.angles.add(new_obj)
                new_obj.atom1.angles.add(new_obj)
                new_obj.atom2.angles.add(new_obj)
                new_obj.atom3.angles.add(new_obj)

        # add dihedral to compound
        elif isinstance(new_obj, Dihedral):
            # don't add A-B-C-D if D-C-B-A is already added
            if not new_obj in self.dihedrals and not new_obj.cloneImage() in self.dihedrals:
                self.dihedrals.add(new_obj)
                new_obj.atom1.dihedrals.add(new_obj)
                new_obj.atom2.dihedrals.add(new_obj)
                new_obj.atom3.dihedrals.add(new_obj)
                new_obj.atom4.dihedrals.add(new_obj)

        # support batch add
        elif isinstance(new_obj, (list, tuple)):
            for elem in new_obj:
                self.add(elem)

        else:
            raise Exception("can't add unknown type " + str(new_obj))


    @staticmethod
    def createEquivalenceTransform(equiv):
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
                for atom0 in pair[0].atoms():
                    self_points = vstack([self_points, atom0.pos])
                for atom1 in pair[1].atoms():
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
        arr = np.fromiter(chain.from_iterable(atom.pos for atom in self.atoms()), dtype=np.float64)
        arrnx3 = arr.reshape((-1,3))

        arrnx3 = T.apply(arrnx3)
        arr = arrnx3.reshape((-1))

        # write back new coordinates into atoms
        i = 0
        for atom in self.atoms():
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

        for part in self.parts:
            # add local atoms
            if isinstance(part, Atom):
                yield part
                # add atoms in sub-components recursively
            if isinstance(part, Compound):
                for subatom in part.atoms():
                    yield subatom


    # def component(self, component_path):
    #     """
    #     Find a component by label
    #     :param component_path: the label of the component (may be hierarchical)
    #     :return: the component (if found), None otherwise
    #     """
    #     dot_pos = component_path.find('.')
    #     if dot_pos > -1:
    #         subcomponent_path = component_path[:dot_pos]
    #         subpath = component_path[dot_pos + 1:]
    #         if subcomponent_path in self._components.keys():
    #             return self._components[subcomponent_path].component(subpath)
    #         else:
    #             return None
    #     else:
    #         if component_path in self._components.keys():
    #             return self._components[component_path]
    #         else:
    #             return None

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

        for a in self.atoms():
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

    def rename(self, array_of_pairs):
        for pair in array_of_pairs:
            for atom in self.getAtomsByKind(pair[0]):
                atom.kind = pair[1]

    # def plot(self, verbose=False, labels=True):
    #     """
    #     Plot atoms in 3d
    #     :param verbose: if True, atom types will be plotted
    #     :param labels: if True, labels will be plotted
    #     """
    #
    #     x = []
    #     y = []
    #     z = []
    #     r = []
    #     g = []
    #     b = []
    #     rgb = []
    #
    #     # sort atoms by type
    #     d = dict()
    #
    #     for (label, atom) in self.atoms():
    #         if atom.kind != 'G' or verbose:
    #             if not atom.kind in d.keys():
    #                 d[atom.kind] = [atom]
    #             else:
    #                 d[atom.kind].append(atom);
    #
    #     for (kind,atomList) in d.items():
    #         x = []
    #         y = []
    #         z = []
    #         r = []
    #         for atom in atomList:
    #             x.append(atom.pos[0])
    #             y.append(atom.pos[1])
    #             z.append(atom.pos[2])
    #             r.append(atom.vdw_radius)
    #
    #         fig = mlab.points3d(x,y,z,r,color=atomList[0].colorRGB, scale_factor=1, scale_mode='scalar')
    #         #fig.glyph.glyph.clamping = False
    #     mlab.show()
    #
    # def plot2(self, verbose=False, labels=False):
    #     """
    #     Plot atoms in 3d
    #     :param verbose: if True, atom types will be plotted
    #     :param labels: if True, labels will be plotted
    #     """
    #     fig = pyplot.figure()
    #     ax = fig.add_subplot(111, projection='3d', aspect='equal')
    #     coord_min = inf
    #     coord_max = -inf
    #     for (label, atom) in self.atoms():
    #         if atom.kind != 'G' or verbose:
    #             # print atom
    #             if labels:
    #                 atom.plot(ax, str(atom))
    #             else:
    #                 atom.plot(ax, None)
    #
    #     ax.set_xlabel('X')
    #     ax.set_ylabel('Y')
    #     ax.set_zlabel('Z')
    #     # ax.set_title(self.label())
    #
    #     pyplot.show()

    # def plot3(self, verbose=False, labels=False):
    #
    #     import networkx as nx
    #     import matplotlib.pyplot as plt
    #
    #     g = nx.Graph()
    #     g.add_node("root")
    #     self.plot3_worker(g, "root")
    #
    #     print "g=" + str(g)
    #
    #     # draw graph
    #     pos = nx.shell_layout(g)
    #     nx.draw(g, pos)
    #
    #     # show graph
    #     plt.show()
    #
    #
    # def plot3_worker(self, g, parent):
    #     print self.__class__
    #
    #     for (label, c) in self._components.iteritems():
    #         mylabel = label+str(c.__hash__())
    #         g.add_node(mylabel)
    #         g.add_edge(parent, mylabel)
    #         print parent + " " + label
    #         if isinstance(c, Compound):
    #             c.plot3_worker(g, mylabel)
    #


    # def __copy__(self):
    #     cls = self.__class__
    #     newone = cls.__new__(cls)
    #
    #     newone.__dict__.update(self.__dict__)
    #     return newone
    #
    # def __deepcopy__(self, memo):
    #     cls = self.__class__
    #     newone = cls.__new__(cls)
    #
    #     memo[id(self)] = newone
    #
    #     newone.parts = set()
    #     newone.references = dict()
    #
    #     for part in self.parts.iteritems():
    #         child_copy = deepcopy(part, memo)
    #         newone.add(child_copy, label)
    #
    #         if label in self._aliases.values():
    #              alias_label = self._aliases.keys()[self._aliases.values().index(label)]
    #              newone.addAlias(alias_label, label)
    #
    #     # copy over the rest of the stuff, that's not yet there in the new one
    #     for k, v in self.__dict__.items():
    #         if not hasattr(newone, k):
    #             setattr(newone, k, deepcopy(v, memo))
    #
    #     return newone




    def getAtomsByKind(self, kind):
        return ifilter(lambda atom: (atom.kind == kind), self.atoms())

    def getBondsByAtomKind(self, kind1, kind2):
        return ifilter(lambda bond: bond.hasAtomKinds(kind1, kind2), self.bonds)

    def getAnglesByAtomKind(self, kind1, kind2, kind3):
        return ifilter(lambda angle: angle.hasAtomKinds(kind1, kind2, kind3), self.angles)

    def getDihedralsByAtomKind(self, kind1, kind2, kind3, kind4):
        return ifilter(lambda dihedral: dihedral.hasAtomKinds(kind1, kind2, kind3, kind4), self.dihedrals)

    def initAtomKdTree(self):
        self.atomsList = list(self.atoms())
        self.atomKdtree = PeriodicCKDTree([atom.pos for atom in self.atomsList], bounds=self.bounds)

    def getAtomsInRange(self, point, radius, maxItems=50):
        # create kdtree if it's not yet there
        if not hasattr(self, 'atomKdtree'):
            self.initAtomKdTree()

        distances, indices = self.atomKdtree.query(point, maxItems)
        # indices = self.kdtree.query(point, maxAtoms)

        neighbors = []
        for index, distance in zip(indices, distances):
            if distance <= radius:
                neighbors.append(self.atomsList[index])
            else:
                break

        return neighbors

    def initBondKdTree(self):
        self.bondsList = list(self.bonds)
        self.BondKdtree = PeriodicCKDTree([bond.com() for bond in self.bondsList], bounds=self.bounds)

    def getBondsInRange(self, point, radius, maxItems=50):
        # create kdtree if it's not yet there
        if not hasattr(self, 'bondKdtree'):
            self.initBondKdTree()

        distances, indices = self.BondKdtree.query(point, maxItems)
        # indices = self.kdtree.query(point, maxAtoms)

        neighbors = []
        for index, distance in zip(indices, distances):
            if distance <= radius:
                neighbors.append(self.bondsList[index])
            else:
                break

        return neighbors

    def initAngleKdTree(self):
        self.anglesList = list(self.angles)
        self.AngleKdtree = PeriodicCKDTree([angle.atom2.pos for angle in self.anglesList], bounds=self.bounds)

    def getAnglesInRange(self, point, radius, maxItems=50):
        # create kdtree if it's not yet there
        if not hasattr(self, 'angleKdtree'):
            self.initAngleKdTree()

        distances, indices = self.AngleKdtree.query(point, maxItems)
        # indices = self.kdtree.query(point, maxAtoms)

        neighbors = []
        for index, distance in zip(indices, distances):
            if distance <= radius:
                neighbors.append(self.anglesList[index])
            else:
                break

        return neighbors

    def findMissingAngleKinds(self):
        missing = set()
        for abc in self.findMissingAngles():
            missing.add((abc.atom1.kind, abc.atom2.kind, abc.atom3.kind))
        return missing

    def findMissingAngles(self):

        missing = set()

        for ab in self.bonds:
            assert(isinstance(ab,Bond))
            type_A = ab.atom1.__class__
            type_B = ab.atom2.__class__

            for bc in ab.atom2.bonds:
                if ab==bc:
                    continue
                abc = Angle.createFromBonds(ab,bc)
                abcImage = abc.cloneImage();
                if not abc in self.angles and not abc.cloneImage() in self.angles and not abcImage in missing:
                    if abc.atom1.kind < abc.atom3.kind:
                        missing.add(abc)
                    else:
                        missing.add(abcImage)

            for bc in ab.atom1.bonds:
                if ab==bc:
                    continue
                abc = Angle.createFromBonds(ab,bc)
                abcImage = abc.cloneImage()
                if not abc in self.angles and not abc.cloneImage() in self.angles and not abcImage in missing:
                    if abc.atom1.kind < abc.atom3.kind:
                        missing.add(abc)
                    else:
                        missing.add(abcImage)

        return missing

    def min_periodic_distance(self, pos1, pos2):
        p1 = np.array(pos1)
        p2 = np.array(pos2)

        # print p1
        # print p2
        p2_img = np.array([[p2[0]+( 0)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           \
                           [p2[0]+( 0)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+( 0)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           \
                           [p2[0]+( 0)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           [p2[0]+( 0)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           \
                           [p2[0]+( 0)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+( 0)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           [p2[0]+( 0)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+( 0)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+( 0)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+( 0)*self.bounds[2]], \
                           \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+( 1)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+( 1)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+( 1)*self.bounds[2]], \
                           [p2[0]+(-1)*self.bounds[0], p2[1]+(-1)*self.bounds[1], p2[2]+(-1)*self.bounds[2]], \
                           ])


        # print p2_img

        # print p2_img-p1
        dists = norm(p2_img-p1, axis=1)
        # print dists

        return dists.min()


# if __name__ == "__main__":
#     c = Compound(bounds=[10,10,10])
#     print c.min_periodic_distance((1,1,1),(9,9,9))