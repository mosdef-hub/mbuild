from collections import OrderedDict, defaultdict
from copy import deepcopy
from orderedset import OrderedSet
from itertools import *

import numpy as np
from numpy.linalg import norm

from atom import *
from angle import Angle
from bond import Bond
from dihedral import Dihedral
from periodic_kdtree import PeriodicCKDTree
from coordinate_transform import *
from warnings import warn

class Compound(object):

    def __init__(self, ctx={}, kind=None, periodicity = None):
        # set the context (optional)
        self.ctx = ctx

        if not periodicity:
            periodicity = np.array([0.0, 0.0, 0.0])

        self.periodicity = periodicity
        # set the kind (defaults to classname)
        if kind:
            self.kind = kind
        else:
            self.kind = self.__class__.__name__

        self.parent = None
        self.parts = OrderedSet() # contains children (compounds or atoms)
        self.references = OrderedDict()  # label to compound/atom mappings -- need not be in parts

        self.bonds = OrderedSet()
        self.angles = OrderedSet()
        self.dihedrals = OrderedSet()

        self.treeDistancePenalty = 1

    def add(self, new_obj, label=None, containment=True, replace=False, inherit_periodicity=True):
        if containment:
            # add atom as a part
            if isinstance(new_obj, Atom):
                self.parts.add(new_obj)
                assert(not new_obj.parent)
                new_obj.parent = self

            # add compound as a part
            elif isinstance(new_obj, Compound):
                self.parts.add(new_obj)
                assert(not new_obj.parent)
                new_obj.parent = self

            # add bond to compound
            elif isinstance(new_obj, Bond):
                # don't add B-A if A-B is already added
                # print "self.bonds before adding: " + str(list(self.bonds))
                # print "adding bond "+str(new_obj)
                if not new_obj in self.bonds and not new_obj.cloneImage() in self.bonds:
                    if containment:
                        self.bonds.add(new_obj)
                    new_obj.atom1.bonds.add(new_obj)
                    new_obj.atom2.bonds.add(new_obj)

            # add angle to compound
            elif isinstance(new_obj, Angle):
                # don't add A-B-C if C-B-A is already added
                if not new_obj in self.angles and not new_obj.cloneImage() in self.angles:
                    if containment:
                        self.angles.add(new_obj)
                    new_obj.atom1.angles.add(new_obj)
                    new_obj.atom2.angles.add(new_obj)
                    new_obj.atom3.angles.add(new_obj)

            # add dihedral to compound
            elif isinstance(new_obj, Dihedral):
                # don't add A-B-C-D if D-C-B-A is already added
                if not new_obj in self.dihedrals and not new_obj.cloneImage() in self.dihedrals:
                    if containment:
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

        # # autogenerate a label for if a reference-only object does not have one
        # if not containment and label is None:
        #     label = new_obj.kind  + "__" + str(id(new_obj))

        if isinstance(new_obj, Compound) or isinstance(new_obj, Atom):
            if label is None:
                label = new_obj.kind  + "__" + str(id(new_obj))

            # add new_obj to references
            if label is not None:
                # support label with counter
                if label.endswith("#"):
                    i = 0
                    while label.replace("#", str(i)) in self.references.keys():
                        i += 1
                    label = label.replace("#", str(i))

                if not replace and label in self.references.keys():
                    raise Exception("label " + label + " already exists in " + str(self))
                else:
                    self.references[label] = new_obj
                    setattr(self, label, new_obj)

            if inherit_periodicity:
                if isinstance(new_obj, Compound):
                    if np.any(new_obj.periodicity):
                        if np.any(self.periodicity):
                            warn("Overriding periodicity of component " + str(self))
                        self.periodicity = new_obj.periodicity

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
            atom.pos = np.array([arr[i], arr[i+1], arr[i+2]])
            i=i+3

        return self


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


    def atomKinds(self):
        kinds = set()
        for atom in self.atoms():
            kinds.add(atom.kind)
        return kinds

    def hasAtomKind(self, kind):
        return kind in self.atomKinds()

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

        min = np.array([minx, miny, minz])
        max = np.array([maxx, maxy, maxz])


        return min, max, max-min

    def rename(self, array_of_pairs):
        for pair in array_of_pairs:
            for atom in self.getAtomListByKind(pair[0]):
                atom.kind = pair[1]

    # def boundingbox_diameter(self, excludeG=True):
    #     (minx, miny, minz), (maxx, maxy, maxz) = self.boundingbox(excludeG=excludeG)
    #     return max([maxx-minx, maxy-miny, maxz-minz])

    def __copy__(self):
        cls = self.__class__
        newone = cls.__new__(cls)

        newone.__dict__.update(self.__dict__)

        newone.parent = None

        return newone

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)

        memo[id(self)] = newone

        for k, v in self.__dict__.items():
            if k is "parent":
                if id(self.parent) in memo:
                    newone.parent = memo[id(self.parent)]
                else:
                    newone.parent = None
            else:
                setattr(newone, k, deepcopy(v, memo))

        return newone

    def initAtomsByKind(self, kind='*'):
        # remember the hash of the parts dict at time of generating the atomListsByKind dict
        self.atomListsByKind_hash = self.parts.__hash__
        self.atomListsByKind = OrderedDict()

        # print "intiializing atoms by kind dict"

        self.atomListsByKind['*'] = []

        for atom in self.atoms():

            self.atomListsByKind['*'].append(atom)

            if atom.kind not in self.atomListsByKind:
                self.atomListsByKind[atom.kind] = [atom]
            else:
                self.atomListsByKind[atom.kind].append(atom)

    def hasAtomListByKind(self, kind='*'):
        if not hasattr(self, 'atomListsByKind') or self.parts.__hash__ != self.atomListsByKind_hash:
            # print "nonexistent of outdated atomsListByKind"
            return False
        else:
            return True


    def getAtomListByKind(self, kind='*'):
        # this is slow...
        # return ifilter(lambda atom: (atom.kind == kind), self.atoms())

        # use some precomputed data structures instead (memory vs. time trade-off)
        if not self.hasAtomListByKind(kind):
            self.initAtomsByKind(kind)

        if kind in self.atomListsByKind:
            return self.atomListsByKind[kind]
        else:
            return []

    def resetAtomListByKind(self, kind='*'):
        if hasattr(self, 'atomListsByKind'):
            del self.atomListsByKind
            del self.atomListsByKind_hash

    def getAtomsByBondType(self, bond_type):
        return ifilter(lambda atom: (atom.bond_type == bond_type), self.atoms())

    def initBondsByAtomKind(self):
        # remember the hash of the bonds dict at time of generating the bondsByAtomKind dict
        self.bondsByAtomKind_hash = self.bonds.__hash__
        self.bondsByAtomKind = dict()

        # print "intiializing bonds by atom kind dict"

        for bond in self.bonds:
            if bond.atom1.kind < bond.atom2.kind:
                pair = (bond.atom1.kind,bond.atom2.kind)
            else:
                pair = (bond.atom2.kind,bond.atom1.kind)

            if pair not in self.bondsByAtomKind:
                self.bondsByAtomKind[pair] = [bond]
            else:
                self.bondsByAtomKind[pair].append(bond)

    def hasBondsByAtomKind(self):
        # check if helper dict doesn't exist, or if the bonds have changed since last computation
        if not hasattr(self, 'bondsByAtomKind_hash') or self.bonds.__hash__ != self.bondsByAtomKind_hash:
            return False
        else:
            return True

    def getBondsByAtomKind(self, kind1, kind2):
        # this runs slowly..
        # return ifilter(lambda bond: bond.hasAtomKinds(kind1, kind2), self.bonds)

        # this should speed it up...
        if not self.hasBondsByAtomKind(): # precompute a data structure if needed
            self.initBondsByAtomKind()

        if kind1 < kind2:
            pair = (kind1, kind2)
        else:
            pair = (kind2, kind1)

        if pair in self.bondsByAtomKind:
            return self.bondsByAtomKind[pair]
        else:
            return []

    def initAnglesByAtomKind(self):
        # remember the hash of the bonds dict at time of generating the anglesByAtomKind dict
        self.anglesByAtomKind_hash = self.angles.__hash__
        self.anglesByAtomKind = dict()

        # print "intializing angles by atom kind dict"

        for angle in self.angles:
            if angle.atom1.kind < angle.atom3.kind:
                triplet = (angle.atom1.kind,angle.atom2.kind,angle.atom3.kind)
            else:
                triplet = (angle.atom3.kind,angle.atom2.kind,angle.atom1.kind)

            if triplet not in self.anglesByAtomKind:
                self.anglesByAtomKind[triplet] = [angle]
            else:
                self.anglesByAtomKind[triplet].append(angle)


    def hasAnglesByAtomKind(self):
        # check if helper dict doesn't exist, or if the bonds have changed since last computation
        if not hasattr(self, 'anglesByAtomKind_hash') or self.bonds.__hash__ != self.anglesByAtomKind_hash:
            return False
        else:
            return True

    def getAnglesByAtomKind(self, kind1, kind2, kind3):

        # this is slow...
        #return ifilter(lambda angle: angle.hasAtomKinds(kind1, kind2, kind3), self.angles)

        # this should speed it up...
        if not self.hasAnglesByAtomKind(): # precompute a data structure if needed
            self.initAnglesByAtomKind()

        if kind1 < kind3:
            triplet = (kind1, kind2, kind3)
        else:
            triplet = (kind3, kind2, kind1)

        if triplet in self.anglesByAtomKind:
            return self.anglesByAtomKind[triplet]
        else:
            return []


    def getDihedralsByAtomKind(self, kind1, kind2, kind3, kind4):
        return ifilter(lambda dihedral: dihedral.hasAtomKinds(kind1, kind2, kind3, kind4), self.dihedrals)

    def initAtomKdTree(self, kind='*'):
        # check if atomKdTrees dict exists and is up-to-date
        if not hasattr(self, 'atomKdTrees') or self.parts.__hash__ != self.atomKdTrees_hash:
            # remember the hash of the bonds dict at time of generating the bondsByAtomKind dict
            self.atomKdTrees_hash = self.parts.__hash__
            self.atomKdTrees = dict()
            # print "intiializing atomKdTrees dict"

        self.atomKdTrees[kind] = PeriodicCKDTree([atom.pos for atom in self.getAtomListByKind(kind)], bounds=self.periodicity)

    def hasAtomKdTree(self, kind='*'):
        if not hasattr(self, 'atomKdTrees') or self.parts.__hash__ != self.atomKdTrees_hash:
            return False

        if kind in self.atomKdTrees:
            return True

        return False

    def getAtomKdTree(self, kind='*'):
        return self.atomKdTrees[kind]

    def getAtomsInRange(self, point, radius, maxItems=10, kind='*'):

        # create kdtree if it's not yet there
        if not self.hasAtomKdTree(kind):
            self.initAtomKdTree(kind)

        distances, indices = self.getAtomKdTree(kind).query(point, maxItems)
        # indices = self.kdtree.query(point, maxAtoms)

        neighbors = []
        for index, distance in zip(indices, distances):
            if distance <= radius:
                neighbors.append(self.getAtomListByKind(kind)[index])
            else:
                break

        return neighbors

    def initBondKdTree(self):
        self.bondsList = list(self.bonds)
        self.bond_kdtree = PeriodicCKDTree([bond.com() for bond in self.bondsList], bounds=self.periodicity)

    def getBondsInRangeKdTree(self, point, radius, maxItems=50):
        # create kdtree if it's not yet there
        if not hasattr(self, 'bond_kdtree'):
            self.initBondKdTree()

        distances, indices = self.bond_kdtree.query(point, maxItems)
        # indices = self.kdtree.query(point, maxAtoms)

        neighbors = []
        for index, distance in zip(indices, distances):
            if distance <= radius:
                neighbors.append(self.bondsList[index])
            else:
                break

        return neighbors

    def getBondsInRange(self, bond, depth=1):
        neighbors = [] # bonds
        for b in bond.atom1.bonds:
            if b is not bond:
                neighbors.append(b)

        for b in bond.atom2.bonds:
            if b is not bond:
                neighbors.append(b)

        return neighbors

    def initAngleKdTree(self):
        self.anglesList = list(self.angles)
        self.AngleKdtree = PeriodicCKDTree([angle.atom2.pos for angle in self.anglesList], bounds=self.periodicity)

    def getAnglesInRange(self, angle, depth=1):
        neighbors = [] # angles
        for a in angle.atom1.angles:
            if a is not angle:
                neighbors.append(a)
        for a in angle.atom2.angles:
            if a is not angle:
                neighbors.append(a)
        for a in angle.atom3.angles:
            if a is not angle:
                neighbors.append(a)

        return neighbors


    def getAnglesInRangeKdTree(self, point, radius, maxItems=50):
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
                abcImage = abc.cloneImage()
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

    def min_periodic_distance(self, x0, x1):
        """Vectorized distance calculation considering minimum image
        """
        d = np.abs(x0 - x1)
        d = np.where(d > 0.5 * self.periodicity, self.periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

    def unique_types(self, ff):
        """
        """
        a_types = OrderedDict()
        i = 1
        for atom in self.atoms():
            if atom.kind != 'G':
                if atom.kind not in a_types:
                    a_types[atom.kind] = (i, ff.atom_types[atom.kind])
                    i += 1
        
        b_types = OrderedDict()
        i = 1
        for bond in self.bonds:
            if bond.kind not in b_types:
                pair = tuple(bond.kind.split('-'))
                b_types[bond.kind] = (i, ff.bond_types[pair])
                i += 1
    
        ang_types = OrderedDict()
        i = 1
        for ang in self.angles:
            if ang.kind not in ang_types:
                triplet = tuple(ang.kind.split('-'))
                ang_types[ang.kind] = (i, ff.angle_types[triplet])
                i += 1

        dih_types = OrderedDict()
        i = 1
        for dih in self.dihedrals:
            if dih.kind not in dih_types:
                quadruplet = tuple(dih.kind.split('-'))
                dih_types[dih.kind] = (i, ff.dih_types[quadruplet])
                i += 1

        return a_types, b_types, ang_types, dih_types

