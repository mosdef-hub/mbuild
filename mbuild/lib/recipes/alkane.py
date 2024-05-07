"""mBuild recipe for a generic alkane chain."""

import mbuild as mb
from mbuild.lib.moieties import CH2, CH3
from mbuild.lib.molecules import Ethane, Methane


class Alkane(mb.Compound):
    """An alkane which may optionally end with a hydrogen or a Port."""

    def __init__(self, n=3, cap_front=True, cap_end=True):
        """Initialize an Alkane Compound.

        Args
        ----
        n: int, default 3,
            Number of carbon atoms.
        cap_front: bool, default True
            Whether to add a methyl group to the beginning of chain
            ('down' port).
        cap_end: bool, default True,
            Whether to add a methyl group to the end of chain ('up' port).
        """
        if n < 1:
            raise ValueError("n must be 1 or more")
        super(Alkane, self).__init__()
        from mbuild.lib.recipes import Polymer

        # Handle the case of Methane and Ethane separately
        if n < 3:
            if n == 1:
                if cap_front and cap_end:
                    self.add(Methane(), "chain")
                elif cap_front != cap_end:
                    chain = CH3()
                    self.add(chain, "chain")
                    if cap_front:
                        self.add(chain["up"], "down", containment=False)
                    else:
                        self.add(chain["up"], "up", containment=False)
                else:
                    chain = CH2()
                    self.add(chain)
                    self.add(chain["down"], "down", containment=False)
                    self.add(chain["up"], "up", containment=False)
            elif n == 2:
                if cap_front and cap_end:
                    self.add(Ethane(), "chain")
                elif cap_front != cap_end:
                    chain = CH2()
                    self.add(chain, "chain")
                    if cap_front:
                        self.add(CH3(), "methyl_front")
                        mb.force_overlap(
                            move_this=self["chain"],
                            from_positions=self["chain"]["up"],
                            to_positions=self["methyl_front"]["up"],
                        )
                        self.add(chain["down"], "down", containment=False)
                    else:
                        self.add(CH3(), "methyl_end")
                        mb.force_overlap(
                            self["methyl_end"],
                            self["methyl_end"]["up"],
                            self["chain"]["down"],
                        )
                        self.add(chain["up"], "up", containment=False)
                else:
                    chain = Polymer(monomers=[CH2()])
                    chain.build(2, add_hydrogens=False)
                    self.add(chain, "chain")
                    self.add(chain["down"], "down", containment=False)
                    self.add(chain["up"], "up", containment=False)

        # Handle general case of n >= 3
        else:
            end_groups = [None, None]
            # Adjust length of Polmyer for absence of methyl terminations.
            if cap_front:
                n -= 1
                end_groups[0] = CH3()
            if cap_end:
                n -= 1
                end_groups[1] = CH3()

            chain = Polymer(monomers=[CH2()], end_groups=end_groups)
            chain.build(n, add_hydrogens=False)
            self.add(chain, "chain")
            if not cap_front:
                # Hoist port label to Alkane level.
                self.add(self["chain"]["up"], "up", containment=False)
            if not cap_end:
                # Hoist port label to Alkane level.
                self.add(self["chain"]["down"], "down", containment=False)
