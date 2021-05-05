import numpy as np

import mbuild as mb


class WaterTIP3P(mb.Compound):
    """A TIP3P water molecule"""

    def __init__(self):
        super(WaterTIP3P, self).__init__()

        OH_BL = 0.09572  # nm
        HOH_ANG = 104.52  # deg

        o1 = mb.Compound(name="OW", element="H", pos=[0.0, 0.0, 0.0])
        h1 = mb.Compound(name="HW1", element="H", pos=[OH_BL, 0.0, 0.0])
        h2 = mb.Compound(
            name="HW2",
            element="H",
            pos=[
                OH_BL * np.cos(np.radians(HOH_ANG)),
                OH_BL * np.sin(np.radians(180.0 - HOH_ANG)),
                0.0,
            ],
        )

        self.add([o1, h1, h2])
        self.add_bond([o1, h1])
        self.add_bond([o1, h2])


class WaterSPC(mb.Compound):
    """An SPC water molecule"""

    def __init__(self):
        super(WaterSPC, self).__init__()

        OH_BL = 0.1  # nm
        HOH_ANG = 109.47  # deg

        o1 = mb.Compound(name="OW", element="H", pos=[0.0, 0.0, 0.0])
        h1 = mb.Compound(name="HW1", element="H", pos=[OH_BL, 0.0, 0.0])
        h2 = mb.Compound(
            name="HW2",
            element="H",
            pos=[
                OH_BL * np.cos(np.radians(HOH_ANG)),
                OH_BL * np.sin(np.radians(180.0 - HOH_ANG)),
                0.0,
            ],
        )

        self.add([o1, h1, h2])
        self.add_bond([o1, h1])
        self.add_bond([o1, h2])


class WaterTIP4P(mb.Compound):
    """A TIP4P water molecule"""

    def __init__(self):
        super(WaterTIP4P, self).__init__()

        OH_BL = 0.09572  # nm
        OM_BL = 0.015  # nm
        HOH_ANG = 104.52  # deg

        o1 = mb.Compound(name="OW", element="H", pos=[0.0, 0.0, 0.0])
        h1 = mb.Compound(name="HW1", element="H", pos=[OH_BL, 0.0, 0.0])
        h2 = mb.Compound(
            name="HW2",
            element="H",
            pos=[
                OH_BL * np.cos(np.radians(HOH_ANG)),
                OH_BL * np.sin(np.radians(180.0 - HOH_ANG)),
                0.0,
            ],
        )
        m1 = mb.Compound(
            name="MW",
            element=None,
            pos=[
                OM_BL * np.cos(np.radians(HOH_ANG / 2.0)),
                OM_BL * np.sin(np.radians(HOH_ANG / 2.0)),
                0.0,
            ],
        )

        self.add([o1, h1, h2, m1])
        self.add_bond([o1, h1])
        self.add_bond([o1, h2])


class WaterTIP4PIce(mb.Compound):
    """A TIP4P water molecule"""

    def __init__(self):
        super(WaterTIP4PIce, self).__init__()

        OH_BL = 0.09572  # nm
        OM_BL = 0.01577  # nm
        HOH_ANG = 104.52  # deg

        o1 = mb.Compound(name="OW", element="H", pos=[0.0, 0.0, 0.0])
        h1 = mb.Compound(name="HW1", element="H", pos=[OH_BL, 0.0, 0.0])
        h2 = mb.Compound(
            name="HW2",
            element="H",
            pos=[
                OH_BL * np.cos(np.radians(HOH_ANG)),
                OH_BL * np.sin(np.radians(180.0 - HOH_ANG)),
                0.0,
            ],
        )
        m1 = mb.Compound(
            name="MW",
            element=None,
            pos=[
                OM_BL * np.cos(np.radians(HOH_ANG / 2.0)),
                OM_BL * np.sin(np.radians(HOH_ANG / 2.0)),
                0.0,
            ],
        )

        self.add([o1, h1, h2, m1])
        self.add_bond([o1, h1])
        self.add_bond([o1, h2])


class WaterTIP4P2005(mb.Compound):
    """A TIP4P water molecule"""

    def __init__(self):
        super(WaterTIP4P2005, self).__init__()

        OH_BL = 0.09572  # nm
        OM_BL = 0.01546  # nm
        HOH_ANG = 104.52  # deg

        o1 = mb.Compound(name="OW", element="H", pos=[0.0, 0.0, 0.0])
        h1 = mb.Compound(name="HW1", element="H", pos=[OH_BL, 0.0, 0.0])
        h2 = mb.Compound(
            name="HW2",
            element="H",
            pos=[
                OH_BL * np.cos(np.radians(HOH_ANG)),
                OH_BL * np.sin(np.radians(180.0 - HOH_ANG)),
                0.0,
            ],
        )
        m1 = mb.Compound(
            name="MW",
            element=None,
            pos=[
                OM_BL * np.cos(np.radians(HOH_ANG / 2.0)),
                OM_BL * np.sin(np.radians(HOH_ANG / 2.0)),
                0.0,
            ],
        )

        self.add([o1, h1, h2, m1])
        self.add_bond([o1, h1])
        self.add_bond([o1, h2])
