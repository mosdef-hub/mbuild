from __builtin__ import classmethod
from itertools import *

from mbuild.xyz import *
from mbuild.bond import *
from mbuild.angle import *
from scipy.spatial import cKDTree
from scipy import inf

__author__ = 'sallai'





class MoleculeModel(object):
    def __init__(self):
        object.__init__(self)
        self.atoms = set()
        self.bonds = set()
        self.angles = set()

    @classmethod
    def create(cls):
        model = MoleculeModel()
        return model

    def add(self, what):
        if isinstance(what, Atom):
            self.atoms.add(what)
            what.bonds = set()
            what.angles = set()
        elif isinstance(what, Bond):
            self.bonds.add(what)
            what.atom1.bonds.add(what)
            what.atom2.bonds.add(what)
        elif isinstance(what, Angle):
            self.angles.add(what)
            what.atom1.angles.add(what)
            what.atom2.angles.add(what)
            what.atom3.angles.add(what)
        elif isinstance(what, (list, tuple)):
            for elem in what:
                self.add(elem)
        else:
            raise Exception("can't add unknown type " + str(what))

    # def plot(self, verbose=False, labels=True):
    #     fig = pyplot.figure()
    #     ax = fig.add_subplot(111, projection='3d', aspect='equal')
    #     ax.set_xlabel('X')
    #     ax.set_ylabel('Y')
    #     ax.set_zlabel('Z')
    #     # ax.set_title(self.label())
    #
    #     for atom in self.atoms:
    #         if atom.atomType != 'G' or verbose:
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

    def getAtomsByType(self, atomType):
        return ifilter(lambda atom: isinstance(atom, atomType), self.atoms)

    def getBondsByTypes(self, atomType1, atomType2):
        return ifilter(lambda bond: bond.hasTypes(atomType1, atomType2), self.bonds)

    def getAnglesByTypes(self, atomType1, atomType2, atomType3):
        return ifilter(lambda angle: angle.hasTypes(atomType1, atomType2, atomType3), self.angles)


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



    def plot(self, verbose=False, labels=True):

        from mayavi import mlab
        x = []
        y = []
        z = []
        r = []
        c = []

        for atom in self.atoms:
            if atom.atomType != 'G' or verbose:
                # print atom
                x.append(atom.pos[0])
                y.append(atom.pos[1])
                z.append(atom.pos[2])
                r.append(atom.vdw_radius)

        mlab.points3d(x,y,z,r)
        """
        for bond in self.bonds:
            epsilon = 0.3
            pos1 = np.array(bond.atom1.pos)
            pos2 = np.array(bond.atom2.pos)
            v12 = pos2 - pos1 # vector from atom1 to atom2
            d12 = np.linalg.norm(v12) # atom1-atom2 distance
            p1 = pos1 + v12/d12 * epsilon
            p2 = pos1 + v12/d12 * (d12 - epsilon)
            mlab.plot3d([p1[0], p2[0]],[p1[1], p2[1]],[p1[2], p2[2]],
                    tube_radius=0.25, color=bond.color)
        """
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
            p2 = pos1 + v12/d12 * (d12 - epsilon)
            mlab.plot3d([pos1[0], pos2[0]],[pos1[1], pos2[1]],[pos1[2], pos2[2]],
                    tube_radius=0.25, color=angle.color)

            mlab.plot3d([pos2[0], pos3[0]],[pos2[1], pos3[1]],[pos2[2], pos3[2]],
                    tube_radius=0.25, color=angle.color)

        mlab.show()


class RuleEngine(object):
    dmin = 0.9
    dmax = 1.1

    @classmethod
    def create(cls, model):
        re = RuleEngine()
        re.model = model
        return re


    def execute(self):
        self.add_bond(C, H, 1.0, 1.4, "c1xh", (1, 1, 1))
        self.add_bond(C, O, 1.2, 1.6, "c1xc", (0, 1, 0))
        self.add_bond(C, C, 1.2, 1.6, "c1xc", (0, 1, 0))
        #self.add_bond(N, C, 1.2, 1.6, "n1xo", (0, 0, 1))
        self.add_angle(C, C, O, "hxcxc", color=(1, 0.5, 0))
        #self.add_angle(C, C, C, "cxcxc", color=(1, 0, 0))
        #self.add_dihedral(H, C, C, C, "hxcxcxc", color=(0, 0, 0))

    def add_bond(self, type_A, type_B, dmin, dmax, bondType, color=(1,1,1)):
        "Ai-Bj distance is in [dmin, dmax] => add bond A1xB(Ai,Bj) (symmetric)"
        for a1 in self.model.getAtomsByType(type_A):
            nearest = self.model.getAtomsInRange(a1.pos, 2)
            for b1 in nearest:
                if isinstance(b1, type_B) and (dmin <= b1.distance(a1) <= dmax):
                    self.model.add(Bond.create(a1, b1, bondType=bondType, color=color))

    def add_angle(self, type_A, type_B, type_C, angleType, thmin=-Inf, thmax=Inf, color=(1,1,1)):
        """
        """
        for ab1 in self.model.getBondsByTypes(type_A, type_B):
            ab = Bond.orderBond(ab1, type_A, type_B)

            nearest = self.model.getBondsInRange(ab.com(), 3)
            for bc1 in nearest:
                if ab1 == bc1:
                    continue
                bc = Bond.orderBond(bc1, type_B, type_C)
                if not bc:
                    continue
                temp_ang = Angle.createFromBonds(ab, bc, angleType=angleType, color=color)
                if temp_ang:
                    if isinstance(temp_ang.atom2, type_B) and (thmin <= temp_ang.inDegrees() <= thmax):
                        self.model.add(temp_ang)

    def add_dihedral(self, type_A, type_B, type_C, type_D, psimin, psimax, dihedralType, color=(1,1,1)):
        """
        """
        for abc1 in self.model.getAnglesByTypes(type_A, type_B, type_C):
            abc = Angle.orderAngle(abc1, type_A, type_B, type_C)

            nearest = self.model.getAnglesInRange(abc.atom2.pos, 5)
            for bcd1 in nearest:
                if abc1 == bcd1:
                    continue
                bcd = Angle.orderAngle(bcd1, type_B, type_C, type_D)
                if (abc.atom2 == bcd.atom1) and (abc.atom3 == bcd.atom2):
                    self.model.add(Dihedral.createFromAngles(abc, bcd, dihedralType=dihedralType, color=color))

    def c1xc1xc(self):
        for b1, b2 in ifilter(lambda (bond1, bond2): bond1.hasCommonAtomsWith(bond2),
                              combinations(
                                      ifilter(lambda bond: bond.bondType == 'c1xc',
                                              self.model.bonds
                                      ),
                                      2
                              )
        ):
            # create angle
            a = Angle.createFromBonds(b1, b2)
            angle_deg = a.inDegrees()
            if 106 <= angle_deg <= 113: # around 109.5
                a.angleType = 'c1xc1xc108'
                a.color = 'red'
                self.model.add(a)
            elif 118 <= angle_deg <= 122: # around 120
                a.angleType = 'c1xc1xc120'
                a.color = 'blue'
                self.model.add(a)
            else:
                print "bad angle:" + str(angle_deg)


    def c1xc1xc_simple(self):
        # Ci-Cj-Ck and Ci != Ck => add angle C1xC1xC(Ci,Cj,Ck) (symmetric)
        amin = 0
        amax = 1

        # for all c
        for ci in self.model.atoms:
            if isinstance(ci, C):
                # for all c1xc bonds that ci is part of
                for ci1xcj in ci.bonds:
                    if ci1xcj.bondType == 'c1xc':
                        cj = ci1xcj.atom1 if ci1xcj.atom1 != ci else ci1xcj.atom2
                        # for all other c1xc bonds of ci
                        for ci1xck in ci.bonds:
                            if ci1xck.bondType == 'c1xc' and ci1xck.atom1 != cj and ci1xck.atom2 != cj:
                                ck = ci1xck.atom1 if ci1xck.atom1 != ci else ci1xck.atom2
                                # create angle
                                angle_deg = Angle.computeInDegrees(cj, ci, ck)
                                if 106 <= angle_deg <= 113:
                                    a = Angle.create(cj, ci, ck, 'c1xc1xc108')
                                    a.color = 'red'
                                    self.model.add(a)
                                elif 118 <= angle_deg <= 122:
                                    a = Angle.create(cj, ci, ck, 'c1xc1xc120')
                                    a.color = 'blue'
                                    self.model.add(a)
                                else:
                                    print "bad angle:" + str(angle_deg)


                                    # "Ci is bonded to Cj and Ci is bonded to Hk => add angle CC1xH(Cj,Ci,Hk) (symmetric)"


                                    # Hi-Cj-Ck-Hl => add dihedral HCC1xH(Hi,Cj,Ck,Hl) (symmetric)


if __name__ == "__main__":
    # m = Xyz.create("c60.xyz")
    m = Xyz.create("../mpc.xyz")
    mm = MoleculeModel.create()
    mm.add(m.atoms().values())


    r = RuleEngine.create(mm)

    r.execute()
    print 'n_bonds: ' + str(len(r.model.bonds))
    print 'n_angles: ' + str(len(r.model.angles))

    # for a in r.model.angles:
    #     print a.inDegrees()

    #print(mm.getAtomsInRange(mm.atoms.pop().pos,2))
    print "done"
    r.model.plot(labels=False)


    # m.plot(labels=False)


