import numpy as np

from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import translate


class Ester(Compound):
    """A ester group -C(=O)O-. """
    def __init__(self):
        super(Ester, self).__init__(self)

        self.append_from_file('ester.pdb')
        translate(self, -self.C[0])

        self.add(Port(anchor=self.O[1]), 'up')
        translate(self.up, self.O[1] + np.array([0, 0.07, 0]))

        self.add(Port(anchor=self.C[0]), 'down')
        translate(self.down, np.array([0, -0.07, 0]))

if __name__ == '__main__':
    m = Ester()
    m.visualize(show_ports=True)

