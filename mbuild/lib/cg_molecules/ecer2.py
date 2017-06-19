import mbuild as mb
import numpy as np

class ECer2(mb.Compound):
    def __init__(self):
        """Returns a CG eCER2 with the head-to-tail vector pointing in -z.
        """
        super(ECer2, self).__init__()
        mb.load('ecer2.hoomdxml', compound=self, relative_to_module=self.__module__)
        self.periodicity = [0, 0, 0]
        xx = list(self.particles())
        mb.coordinate_transform.z_axis_transform(self,
                new_origin=xx[6], point_on_z_axis=xx[11])
        self.rotate(np.pi, [1, 0, 0])
