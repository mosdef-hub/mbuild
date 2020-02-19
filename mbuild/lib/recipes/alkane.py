import mbuild as mb

from mbuild.lib.moieties import CH2
from mbuild.lib.atoms import H


class Alkane(mb.Compound):
    """An alkane which may optionally end with a hydrogen or a Port."""
    def __init__(self, n=3, cap_front=True, cap_end=True):
        """Initialize an Alkane Compound.

        Args:
            n: Number of carbon atoms.
            cap_front: Add a hydrogen to the beginning of chain ('down' port).
            cap_end: Add a hydrogen to the end of chain ('up' port).
        """
        if n < 1:
            raise ValueError('n must be 1 or more')
        super(Alkane, self).__init__()

        # Adjust length of Polmyer for absence of methyl terminations.
        chain = mb.lib.recipes.Polymer(CH2(), n=n, port_labels=('up', 'down'))
        self.add(chain, 'chain')

        if cap_front:
            self.add(H(), "hydrogen_front")
            mb.force_overlap(move_this=self['chain'],
                             from_positions=self['chain']['up'],
                             to_positions=self['hydrogen_front']['up'])
        else:
            # Hoist port label to Alkane level.
            self.add(chain['up'], 'up', containment=False)

        if cap_end:
            self.add(H(), 'hydrogen_end')
            mb.force_overlap(self['hydrogen_end'], self['hydrogen_end']['up'], self['chain']['down'])
        else:
            # Hoist port label to Alkane level.
            self.add(chain['down'], 'down', containment=False)
