from mbuild.compound import Compound
from mbuild.testing.tools import get_fn


class ECeramideNS(Compound):
    """ """
    def __init__(self):
        super(ECeramideNS, self).__init__()
        self.append_from_file(get_fn('ecer2.hoomdxml'))  # this is a CG model
        #self.append_from_file(get_fn('e-ceramide-ns.pdb'))

if __name__ == "__main__":
    ecerns = ECeramideNS()

    from mbuild.plot import Plot
    Plot(ecerns, verbose=True, atoms=True, bonds=True).show()


