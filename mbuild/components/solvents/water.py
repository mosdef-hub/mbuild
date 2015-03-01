import numpy as np

from mbuild.compound import Compound
from mbuild.tools.tiled_compound import TiledCompound


class Water(Compound):
    """An SPC water box."""
    def __init__(self):
        super(Water, self).__init__()

        self.append_from_file('spc216.pdb',
                                relative_to_module=self.__module__)
        self.periodicity = np.array([1.0, 1.0, 1.0])

if __name__ == "__main__":
    s = Water()
    s.visualize(show_ports=True)
