from mbuild.prototype import Prototype
from rules import Rules
__author__ = 'sallai'

class ForceFieldRules(Rules):

    def __init__(self, ff):
        self.force_field = ff

    def apply(self, compound):
        # set atom prototypes
        for atom in compound.atoms():
            if atom.kind in self.force_field.atom_types:
                atomType = self.force_field.atom_types[atom.kind]

                Prototype(atom.kind,
                        bond_type=atomType.alias,
                        atomic_num=atomType.atomicNumber,
                        mass=atomType.mass,
                        charge=atomType.charge,
                        sigma=atomType.sigma,
                        epsilon=atomType.epsilon)
            else:
                print "Atom kind " + str(atom.kind) + " not found in force field"

        # add bonds
        for (type_A, type_B), bondType in self.force_field.bond_types.iteritems():
            if compound.hasAtomKind(type_A) and compound.hasAtomKind(type_B):
                rmin = bondType.r * 0.8
                rmax = bondType.r * 1.2
                print "pair "+type_A+" "+type_B
                self.add_bond(compound, type_A, type_B, rmin, rmax, "{0}-{1}".format(type_A,type_B))