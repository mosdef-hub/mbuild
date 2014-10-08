from mbuild.compound import Compound
from mbuild.examples.ethane.methyl import Methyl
from mbuild.coordinate_transform import equivalence_transform


class Ethane(Compound):
    """ """

    def __init__(self):
        super(Ethane, self).__init__(kind='Ethane')
        self.add(Methyl(), "m1")
        self.add(Methyl(), "m2")

        equivalence_transform(self.m1, self.m1.up, self.m2.down)


def main():
    ethane = Ethane()
    
    ethane = ethane.to_trajectory()
    ethane.top.find_forcefield_terms()
    ethane.save('ethane.hoomdxml')

if __name__ == "__main__":
    main()
