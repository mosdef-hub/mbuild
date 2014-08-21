from mbuild.bond import Bond

__author__ = 'sallai'


class Rules(object):

    def apply(self, compound):
        raise Exception("RuleEngine must be subclassed, with the execute method implemented in the subclass")

    def add_bond(self, compound, type_A, type_B, dmin, dmax, kind):
        "Ai-Bj distance is in [dmin, dmax] => add bond A1xB(Ai,Bj) (symmetric)"
        print "+ "+type_A+" "+type_B

        for a1 in compound.atom_list_by_kind(type_A):
            nearest = compound.getAtomsInRange(a1.pos, dmax)
            print nearest
            for b1 in nearest:
                print "- "+a1.kind+" "+b1.kind

                if (b1.kind==type_B):
                    print "found pair "+type_A+" "+type_B+" distance="+str(compound.min_periodic_distance(b1.pos, a1.pos))+" min="+str(dmin)+" max="+str(dmax)

                if (b1.kind==type_B) and (dmin <= compound.min_periodic_distance(b1.pos, a1.pos) <= dmax):

                    compound.add(Bond(a1, b1, kind=kind))

    # def add_angle(self, type_A, type_B, type_C, kind, thmin=-Inf, thmax=Inf):
    #     """
    #     """
    #
    #     for ab1 in self.compound.getBondsByAtomKind(type_A, type_B):
    #         ab = ab1.cloneWithOrder(type_A, type_B)
    #         nearest = self.compound.getBondsInRange(ab.com(), 10)
    #         for bc1 in nearest:
    #             if ab1 == bc1:
    #                 continue
    #             if not bc1.hasAtomKinds(type_B, type_C):
    #                 continue
    #
    #             bc = bc1.cloneWithOrder(type_B, type_C)
    #
    #             temp_ang = Angle.createFromBonds(ab, bc, kind=kind)
    #             if temp_ang:
    #                 if (temp_ang.atom2.kind==type_B) and (thmin <= temp_ang.inDegrees() <= thmax):
    #                     self.compound.add(temp_ang)
    #
    # def add_dihedral(self, type_A, type_B, type_C, type_D, dihedralKind):
    #     """
    #     """
    #
    #     for abc1 in self.compound.getAnglesByAtomKind(type_A, type_B, type_C):
    #         abc = abc1.cloneWithOrder(type_A, type_B, type_C)
    #
    #         nearest = self.compound.getAnglesInRange(abc.atom2.pos, 5)
    #         for bcd1 in nearest:
    #             if abc1 == bcd1:
    #                 continue
    #
    #             if not bcd1.hasAtomKinds(type_B, type_C, type_D):
    #                 continue
    #
    #             bcd = bcd1.cloneWithOrder(type_B, type_C, type_D)
    #
    #             temp_dhdr = Dihedral.createFromAngles(abc, bcd, kind=dihedralKind)
    #             if temp_dhdr:
    #                 self.compound.add(temp_dhdr)
