import numpy as np

import mbuild as mb


class Water(mb.Compound):
    """An SPC water box."""
    def __init__(self):
        super(Water, self).__init__()

        mb.load('spc216.pdb', compound=self, relative_to_module=self.__module__)
        self.periodicity = np.array([1.0, 1.0, 1.0])

if __name__ == "__main__":
    wat = Water()
    wat.visualize(show_ports=True)
