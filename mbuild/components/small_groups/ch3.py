from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.coordinate_transform import translate


class Ch3(Compound):
    """A methyl group. """
    def __init__(self):
        super(Ch3, self).__init__(self)

        self.append_from_file('ch3.pdb', relative_to_module=self.__module__)
        translate(self, -self.C[0])  # Move carbon to origin.

        self.add(Port(anchor=self.C[0]), 'up')
        translate(self.up, [0, -0.07, 0])

if __name__ == '__main__':
    m = Ch3()
    m.visualize(show_ports=True)

