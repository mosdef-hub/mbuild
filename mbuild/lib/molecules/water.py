import mbuild as mb
import numpy as np


class WaterTIP3P(mb.Compound):
    """A TIP3P water molecule"""
    def __init__(self):
        super(WaterTIP3P, self).__init__()

        OH_BL = 0.09572  # nm
        HOH_ANG = 104.52  # deg

        O = mb.Compound(name="OW", element="H", pos=[0.0, 0.0, 0.0])
        H1 = mb.Compound(name="HW1", element="H", pos=[0.0, OH_BL, 0])
        H2 = mb.Compound(
            name="HW2", element="H",
            pos=[
                OH_BL * np.cos(np.radians(180.0 - HOH_ANG)),
                OH_BL * np.sin(np.radians(180.0 - HOH_ANG)),
                0.0,
            ]
        )

        self.add([O, H1, H2])
        self.add_bonds([O, H1])
        self.add_bonds([O, H2])


