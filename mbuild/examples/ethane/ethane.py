from mbuild.compound import Compound
from mbuild.components.small_groups.ch3 import Ch3
from mbuild.coordinate_transform import equivalence_transform


class Ethane(Compound):
    """An ethane molecule. """
    def __init__(self):
        """Connect two methyl groups to form an ethane. """
        super(Ethane, self).__init__(kind='Ethane')

        self.add(Ch3(), "methyl1")
        self.add(Ch3(), "methyl2")
        equivalence_transform(self.methyl1, self.methyl1.up, self.methyl2.up)


def main():
    ethane = Ethane()
    return ethane


if __name__ == "__main__":
    ethane = main()
    ethane.visualize(show_ports=True)
    import pdb
    pdb.set_trace()
