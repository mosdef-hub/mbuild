import mbuild as mb
import numpy as np

class FFAC16(mb.Compound):
    def __init__(self):
        """Returns a CG FFA C16 with the head-to-tail vector pointing in -z.
        """
        super(FFAC16, self).__init__()
        mb.load('ffac16.hoomdxml', compound=self, relative_to_module=self.__module__)
        self.periodicity = [0, 0, 0]
        xx = list(self.particles())
        mb.coordinate_transform.z_axis_transform(self,
                new_origin=xx[5], point_on_z_axis=xx[0])
        self.spin(np.pi, [1, 0, 0])
