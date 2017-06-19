import mbuild as mb
import numpy as np

class Water(mb.Compound):
    def __init__(self):
        """Returns a CG FFA C16 with the head-to-tail vector pointing in -z.
        """
        super(Water, self).__init__()
        mb.load('water.hoomdxml', compound=self, relative_to_module=self.__module__)
        self.periodicity = [0, 0, 0]
