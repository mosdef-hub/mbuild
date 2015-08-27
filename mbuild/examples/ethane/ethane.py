import mbuild as mb

from mbuild.lib.moieties import CH3


class Ethane(mb.Compound):
    """An ethane molecule. """
    def __init__(self):
        """Connect two methyl groups to form an ethane. """
        super(Ethane, self).__init__(kind='Ethane')

        self.add(CH3(), "methyl1")
        self.add(CH3(), "methyl2")
        mb.equivalence_transform(self.methyl1, self.methyl1.up, self.methyl2.up)


def main():
    ethane = Ethane()
    return ethane


if __name__ == "__main__":
    ethane = main()
    ethane.visualize(show_ports=True)
