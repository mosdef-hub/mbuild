from __builtin__ import classmethod
from itertools import *

from mbuild.xyz import *
from mbuild.bond import *
from mbuild.angle import *


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

    def getAllAtomsByType(self, atomType):
        return ifilter(lambda atom: isinstance(atom, atomType), self.atoms)


    def plot(self, verbose=False, labels=True):

        from mayavi import mlab
        x = []
        y = []
        z = []

        for atom in self.atoms:
            if atom.atomType != 'G' or verbose:
                # print atom
                x.append(atom.pos[0])
                y.append(atom.pos[1])
                z.append(atom.pos[2])

        s = mlab.points3d(x,y,z)

        for bond in self.bonds:
            epsilon = 0.3
            pos1 = np.array(bond.atom1.pos)
            pos2 = np.array(bond.atom2.pos)
            v12 = pos2 - pos1 # vector from atom1 to atom2
            d12 = np.linalg.norm(v12) # atom1-atom2 distance
            p1 = pos1 + v12/d12 * epsilon
            p2 = pos1 + v12/d12 * (d12 - epsilon)
            mlab.plot3d([p1[0], p2[0]],[p1[1], p2[1]],[p1[2], p2[2]])


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
        # self.c1xc()
        self.c1xh()
        self.c1xc1xc()

    def c1xc(self):
        "Ci-Cj distance is in [dmin, dmax] => add bond C1xC(Ci,Cj) (symmetric)"

        dmin = 1.2
        dmax = 1.7

        for c1, c2 in ifilter(lambda (c1, c2): dmin <= c1.distance(c2) <= dmax,
                              combinations(
                                      ifilter(lambda atom: isinstance(atom, C),
                                              self.model.atoms),
                                      2
                              )
        ):
            self.model.add(Bond.create(c1, c2, bondType="c1xc", color='green'))

    def c1xh(self):
        "Ci-Hj distance is in [dmin, dmax] => add bond C1xH(Ci,Hj) (symmetric)"

        dmin = 1.0
        dmax = 1.5
        # print list(combinations(list(self.model.getAllAtomsByType(C))+ list(self.model.getAllAtomsByType(H)),2))



        for c1, h1 in ifilter(lambda (c1, h1): dmin <= c1.distance(h1) <= dmax and isinstance(c1,C) and isinstance(h1,H),
                                      combinations(list(self.model.getAllAtomsByType(C)) + list(self.model.getAllAtomsByType(H)),
                                      2
                              )
        ):
            self.model.add(Bond.create(c1, h1, bondType="c1xh", color='blue'))


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
    print r.model.bonds

    # for a in r.model.angles:
    #     print a.inDegrees()


    r.model.plot(labels=False)


    # m.plot(labels=False)


