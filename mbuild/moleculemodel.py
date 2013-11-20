from itertools import ifilter
from sets import Set
import numpy
from scipy.spatial.ckdtree import cKDTree
from mbuild.angle import Angle
from mbuild.atom import Atom
from mbuild.bond import Bond
from mbuild.dihedral import Dihedral
import numpy as np
__author__ = 'sallai'


class MoleculeModel(object):
    def __init__(self):
        object.__init__(self)
        self.atoms = set()
        self.bonds = set()
        self.angles = set()
        self.dihedrals = set()


    @classmethod
    def create(cls):
        model = MoleculeModel()
        return model

    def add(self, new_obj):
        if isinstance(new_obj, Atom):
            self.atoms.add(new_obj)
            new_obj.bonds = set()
            new_obj.angles = set()
            new_obj.dihedrals = set()
        elif isinstance(new_obj, Bond):
            if not new_obj in self.bonds and not new_obj.cloneImage() in self.bonds:
                self.bonds.add(new_obj)
                new_obj.atom1.bonds.add(new_obj)
                new_obj.atom2.bonds.add(new_obj)
        elif isinstance(new_obj, Angle):
            if not new_obj in self.angles and not new_obj.cloneImage() in self.angles:
                self.angles.add(new_obj)
                new_obj.atom1.angles.add(new_obj)
                new_obj.atom2.angles.add(new_obj)
                new_obj.atom3.angles.add(new_obj)
        elif isinstance(new_obj, Dihedral):
            if not new_obj in self.dihedrals and not new_obj.cloneImage() in self.dihedrals:
                self.dihedrals.add(new_obj)
                new_obj.atom1.dihedrals.add(new_obj)
                new_obj.atom2.dihedrals.add(new_obj)
                new_obj.atom3.dihedrals.add(new_obj)
                new_obj.atom4.dihedrals.add(new_obj)
        elif isinstance(new_obj, (list, tuple)):
            for elem in new_obj:
                self.add(elem)
        else:
            raise Exception("can't add unknown type " + str(new_obj))

    # def plot(self, verbose=False, labels=True):
    #     fig = pyplot.figure()
    #     ax = fig.add_subplot(111, projection='3d', aspect='equal')
    #     ax.set_xlabel('X')
    #     ax.set_ylabel('Y')
    #     ax.set_zlabel('Z')
    #     # ax.set_title(self.label())
    #
    #     for atom in self.atoms:
    #         if atom.kind != 'G' or verbose:
    #             # print atom
    #             if labels:
    #                 atom.plot(ax, str(atom))
    #             else:
    #                 atom.plot(ax, None)
    #
    #     for bond in self.bonds:
    #         bond.plot(ax)
    #
    #     for angle in self.angles:
    #         angle.plot(ax)
    #
    #     pyplot.show()

    def getAtomsByKind(self, kind):
        return ifilter(lambda atom: isinstance(atom, kind), self.atoms)

    def getBondsByAtomKind(self, kind1, kind2):
        return ifilter(lambda bond: bond.hasAtomKinds(kind1, kind2), self.bonds)

    def getAnglesByAtomKind(self, kind1, kind2, kind3):
        return ifilter(lambda angle: angle.hasAtomKinds(kind1, kind2, kind3), self.angles)

    def getDihedralsByAtomKind(self, kind1, kind2, kind3, kind4):
        return ifilter(lambda dihedral: dihedral.hasAtomKinds(kind1, kind2, kind3, kind4), self.dihedrals)

    def initAtomKdTree(self):
        self.atomsList = list(self.atoms)
        self.atomKdtree = cKDTree([atom.pos for atom in self.atomsList])

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
        self.BondKdtree = cKDTree([bond.com() for bond in self.bondsList])

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
        self.AngleKdtree = cKDTree([angle.atom2.pos for angle in self.anglesList])

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



    def plot(self, verbose=False, labels=True, atoms=True, bonds=True, angles=True, dihedrals=True):

        from mayavi import mlab

        # display atoms
        if atoms:
            # sort atoms by type
            d = dict()

            for atom in self.atoms:
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
                    r.append(atom.vdw_radius/10)

                fig = mlab.points3d(x,y,z,r,color=atomList[0].colorRGB, scale_factor=1, scale_mode='scalar')



        if bonds:
            for bond in self.bonds:
                epsilon = 0.4
                pos1 = np.array(bond.atom1.pos)
                pos2 = np.array(bond.atom2.pos)
                v12 = pos2 - pos1 # vector from atom1 to atom2
                d12 = np.linalg.norm(v12) # atom1-atom2 distance
                p1 = pos1 + v12/d12 * epsilon
                p2 = pos1 + v12/d12 * (d12 - epsilon)
                mlab.plot3d([p1[0], p2[0]],[p1[1], p2[1]],[p1[2], p2[2]],
                        tube_radius=0.25, color=bond.color, opacity=1)

        if angles:
            for angle in self.angles:
                epsilon = 0.3
                pos1 = np.array(angle.atom1.pos)
                pos2 = np.array(angle.atom2.pos)
                pos3 = np.array(angle.atom3.pos)
                v12 = pos2 - pos1 # vector from atom1 to atom2
                v23 = pos3 - pos2
                d12 = np.linalg.norm(v12) # atom1-atom2 distance
                d23 = np.linalg.norm(v23) # atom1-atom2 distance
                p1 = pos1 + v12/d12 * epsilon
                # p2 = pos1 + v12/d12 * (d12 - epsilon)
                p3 = pos3 - v23/d23 * epsilon
                mlab.plot3d([pos2[0], p1[0]],[pos2[1], p1[1]],[pos2[2], p1[2]],
                        tube_radius=0.20, color=angle.color, opacity=.5)

                mlab.plot3d([pos2[0], p3[0]],[pos2[1], p3[1]],[pos2[2], p3[2]],
                        tube_radius=0.20, color=angle.color, opacity=.5)

        if dihedrals:
            for d in self.dihedrals:
                epsilon = 0.2
                pos1 = np.array(d.atom1.pos)
                pos2 = np.array(d.atom2.pos)
                pos3 = np.array(d.atom3.pos)
                pos4 = np.array(d.atom4.pos)
                v12 = pos2 - pos1 # vector from atom1 to atom2
                v34 = pos4 - pos3
                d12 = np.linalg.norm(v12) # atom1-atom2 distance
                d34 = np.linalg.norm(v34) # atom1-atom2 distance
                p1 = pos1 + v12/d12 * epsilon
                p4 = pos3 + v34/d34 * (d34 - epsilon)
                mlab.plot3d([pos2[0], p1[0]],[pos2[1], p1[1]],[pos2[2], p1[2]],
                        tube_radius=0.10, color=d.color)

                mlab.plot3d([pos2[0], pos3[0]],[pos2[1], pos3[1]],[pos2[2], pos3[2]],
                        tube_radius=0.10, color=d.color)

                mlab.plot3d([pos3[0], p4[0]],[pos3[1], p4[1]],[pos3[2], p4[2]],
                        tube_radius=0.10, color=d.color)

        mlab.show()