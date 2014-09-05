import time

from mbuild.coordinate_transform import *
from mbuild.file_formats.xyz import Xyz
from mbuild.compound import Compound
from mbuild.ff.opls_rules import OplsRules
from mbuild.ff.opls_forcefield import OplsForceField

class Buckyball(Compound):

    def __init__(self, ctx={}, alpha=0):
        super(Buckyball, self).__init__(ctx=ctx)

        # Look for xyz file in same directory as this file.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        xyz_path = os.path.join(current_dir, 'c60.xyz')

        # Read xyz file with labels.
        # First C is labeled C_0, second C is C_1 etc.
        buckyball = Xyz(xyz_path, labels=True)
        self.add(buckyball, 'buckyball_xyz')

if __name__ == "__main__":
    print "Generating model..."
    start = time.time()
    m = Buckyball()
    print "Done. ({0:.2f} s)".format(time.time() - start)

    print "Loading and pruning forcefield..."
    start = time.time()
    ff = OplsForceField()
    ff = ff.prune(m)
    print "Done. ({0:.2f} s)".format(time.time() - start)

    print "Generating topology..."
    start = time.time()
    rules = OplsRules(m, ff)
    rules.execute()
    print "Done. ({0:.2f} s)".format(time.time() - start)

    print "\nNumber of atoms: {0}".format(len(m.atom_list_by_kind('*')))
    print "Number of bonds: {0}".format(len(m.bonds))
    print "Number of angles: {0}".format(len(m.angles))
    print "Number of dihedrals: {0}".format(len(m.dihedrals))

    from mbuild.plot import Plot
    Plot(m).show()
