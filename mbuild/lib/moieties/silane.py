import numpy as np

import mbuild as mb


class Silane(mb.Compound):
    """An Si(OH)2 group with two ports. """
    def __init__(self, ):
        super(Silane, self).__init__()
        mb.load('silane.pdb', compound=self, relative_to_module=self.__module__)

        # Transform the coordinate system such that the silicon atom is at the
        # origin and the oxygen atoms are on the x axis.
        mb.x_axis_transform(self, new_origin=self[0], point_on_x_axis=self[1])

        # Add bottom port.
        self.add(mb.Port(anchor=self[0]), 'down')
        self['down'].translate(np.array([0, -.07, 0]))

        # Add top port.
        self.add(mb.Port(anchor=self[0]), 'up')
        self['up'].translate(np.array([0, .07, 0]))

if __name__ == "__main__":
    m = Silane()
    m.save('silane.mol2')
