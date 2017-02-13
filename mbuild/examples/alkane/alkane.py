# -- ==alkane== --
import mbuild as mb

from mbuild.lib.moieties import CH2
from mbuild.lib.moieties import CH3


class Alkane(mb.Compound):
    """An alkane which may optionally end with a hydrogen or a Port."""
    def __init__(self, n=3, cap_front=True, cap_end=True):
        """Initialize an Alkane Compound.

        Args:
            n: Number of carbon atoms.
            cap_front: Add methyl group to beginning of chain ('down' port).
            cap_end: Add methyl group to end of chain ('up' port).
        """
        if n < 2:
            raise ValueError('n must be 1 or more')
        super(Alkane, self).__init__()

        # Adjust length of Polmyer for absence of methyl terminations.
        if not cap_front:
            n += 1
        if not cap_end:
            n += 1
        chain = mb.Polymer(CH2(), n=n-2, port_labels=('up', 'down'))
        self.add(chain, 'chain')

        if cap_front:
            self.add(CH3(), "methyl_front")
            mb.force_overlap(move_this=self['chain'],
                             from_positions=self['chain']['up'],
                             to_positions=self['methyl_front']['up'])
        else:
            # Hoist port label to Alkane level.
            self.add(chain['up'], 'up', containment=False)

        if cap_end:
            self.add(CH3(), 'methyl_end')
            mb.force_overlap(self['methyl_end'], self['methyl_end']['up'], self['chain']['down'])
        else:
            # Hoist port label to Alkane level.
            self.add(chain['down'], 'down', containment=False)

# -- ==alkane== --