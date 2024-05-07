"""An ethane molecule."""

import mbuild as mb
from mbuild.lib.moieties import CH3


class Ethane(mb.Compound):
    """An ethane molecule.

    Connect two methyl groups to form an ethane.
    """

    def __init__(self):
        super(Ethane, self).__init__()

        self.add(CH3(), "methyl1")
        self.add(CH3(), "methyl2")
        mb.force_overlap(
            self["methyl1"], self["methyl1"]["up"], self["methyl2"]["up"]
        )
