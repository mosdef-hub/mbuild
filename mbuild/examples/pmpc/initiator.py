import numpy as np

import mbuild as mb


class Initiator(mb.Compound):
    def __init__(self):
        super(Initiator, self).__init__()

        # Look for data file in same directory as this python module.
        mb.load('initiator.pdb', compound=self, relative_to_module=self.__module__)

        # Transform the coordinate system such that the two carbon atoms
        # that are part of the backbone are on the y axis, C_1 at the origin.
        mb.y_axis_transform(self, new_origin=self[0], point_on_y_axis=self[21])

        # Add bottom port
        self.add(mb.Port(anchor=self[0]), 'down')
        # Place the port.
        self['down'].translate(np.array([0.0, -0.07, 0.0]))

        # Add top port.
        self.add(mb.Port(anchor=self[21]), 'up')
        # Place the port.
        self['up'].translate(np.array([0.0, 0.07, 0.0]))

if __name__ == "__main__":
    ini = Initiator()
    print(ini)
