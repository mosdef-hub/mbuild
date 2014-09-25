from mbuild.compound import Compound
from mbuild.testing.tools import get_fn


class ECerNS(Compound):
    """ """
    def __init__(self):
        Compound.__init__(self)
        self.append_from_file(get_fn('ecer2.hoomdxml'))

if __name__ == "__main__":
    ecerns = ECerNS()

    from mbuild.plot import Plot
    Plot(ecerns, verbose=True, atoms=True, bonds=True).show()


