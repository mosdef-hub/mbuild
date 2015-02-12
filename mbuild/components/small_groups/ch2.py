import numpy as np

from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import translate


class Ch2(Compound):
    """A methylene bridge. """
    def __init__(self):
        super(Ch2, self).__init__(self)

        self.append_from_file('ch2.pdb', relative_to_module=self.__module__)

        self.add(Port(anchor=self.C[0]), 'up')
        translate(self.up, np.array([0, 0.07, 0]))

        self.add(Port(anchor=self.C[0]), 'down')
        translate(self.down, np.array([0, -0.07, 0]))

if __name__ == '__main__':
    m = Ch2()
    m.visualize(show_ports=True)

