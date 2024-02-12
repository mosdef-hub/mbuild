"""Water molecules with geometries from different models."""

import numpy as np

import mbuild as mb


class Water3Site(mb.Compound):
    """A generic 3-site water model."""

    def __init__(self, oh_bond_length, hoh_angle):
        super().__init__()

        o1 = mb.Compound(name="OW", element="O", pos=[0.0, 0.0, 0.0])
        h1 = mb.Compound(
            name="HW1", element="H", pos=[oh_bond_length, 0.0, 0.0]
        )
        h2 = mb.Compound(
            name="HW2",
            element="H",
            pos=[
                oh_bond_length * np.cos(np.radians(hoh_angle)),
                oh_bond_length * np.sin(np.radians(180.0 - hoh_angle)),
                0.0,
            ],
        )

        self.add([o1, h1, h2])
        self.add_bond([o1, h1])
        self.add_bond([o1, h2])


class WaterTIP3P(Water3Site):
    """A TIP3P water molecule.

    Paper: https://doi.org/10.1063/1.445869
    Additional reference: https://lammps.sandia.gov/doc/Howto_tip3p.html
    """

    def __init__(self):
        oh_bond_length = 0.09572  # nm
        hoh_angle = 104.52  # deg
        super().__init__(oh_bond_length, hoh_angle)


class WaterSPC(Water3Site):
    """An SPC water molecule.

    Paper: https://doi.org/10.1021/j100308a038
    Additional reference: https://lammps.sandia.gov/doc/Howto_spc.html
    """

    def __init__(self):
        oh_bond_length = 0.1  # nm
        hoh_angle = 109.47  # deg
        super().__init__(oh_bond_length, hoh_angle)


class Water4Site(mb.Compound):
    """A generic 4-site water model."""

    def __init__(self, oh_bond_length, hoh_angle, om_bond_length):
        super().__init__()

        o1 = mb.Compound(name="OW", element="O", pos=[0.0, 0.0, 0.0])
        h1 = mb.Compound(
            name="HW1", element="H", pos=[oh_bond_length, 0.0, 0.0]
        )
        h2 = mb.Compound(
            name="HW2",
            element="H",
            pos=[
                oh_bond_length * np.cos(np.radians(hoh_angle)),
                oh_bond_length * np.sin(np.radians(180.0 - hoh_angle)),
                0.0,
            ],
        )
        m1 = mb.Compound(
            name="MW",
            element=None,
            pos=[
                om_bond_length * np.cos(np.radians(hoh_angle / 2.0)),
                om_bond_length * np.sin(np.radians(hoh_angle / 2.0)),
                0.0,
            ],
        )

        self.add([o1, h1, h2, m1])
        self.add_bond([o1, h1])
        self.add_bond([o1, h2])


class WaterTIP4P(Water4Site):
    """A TIP4P water molecule.

    Paper: https://doi.org/10.1063/1.445869
    Additional reference: https://lammps.sandia.gov/doc/Howto_tip4p.html
    """

    def __init__(self):
        oh_bond_length = 0.09572  # nm
        om_bond_length = 0.015  # nm
        hoh_angle = 104.52  # deg
        super().__init__(oh_bond_length, hoh_angle, om_bond_length)


class WaterTIP4PIce(Water4Site):
    """A TIP4P/Ice water molecule.

    Paper: https://doi.org/10.1063/1.1931662
    Additional reference: https://lammps.sandia.gov/doc/Howto_tip4p.html
    """

    def __init__(self):
        oh_bond_length = 0.09572  # nm
        om_bond_length = 0.01577  # nm
        hoh_angle = 104.52  # deg
        super().__init__(oh_bond_length, hoh_angle, om_bond_length)


class WaterTIP4P2005(Water4Site):
    """A TIP4P/2005 water molecule.

    Paper: https://doi.org/10.1063/1.2121687
    Additional reference: https://lammps.sandia.gov/doc/Howto_tip4p.html
    """

    def __init__(self):
        oh_bond_length = 0.09572  # nm
        om_bond_length = 0.01546  # nm
        hoh_angle = 104.52  # deg
        super().__init__(oh_bond_length, hoh_angle, om_bond_length)
