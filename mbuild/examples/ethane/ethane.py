from mbuild.compound import Compound
from mbuild.examples.ethane.methyl import Methyl
from mbuild.coordinate_transform import equivalence_transform


class Ethane(Compound):
    """An ethane molecule. """
    def __init__(self):
        """Connect two methyl groups to form an ethane. """
        super(Ethane, self).__init__(kind='Ethane')

        self.add(Methyl(), "m1")
        self.add(Methyl(), "m2")
        equivalence_transform(self.m1, self.m1.up, self.m2.down)


def main():
    ethane = Ethane()
    return ethane


if __name__ == "__main__":
    ethane = main()
    import pdb
    pdb.set_trace()
