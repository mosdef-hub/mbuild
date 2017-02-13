from numpy import pi

import mbuild as mb

from mbuild.lib.moieties import Silane
from mbuild.lib.moieties import CH3
from mbuild.examples.pmpc.mpc import MPC
from mbuild.examples.pmpc.initiator import Initiator


class Brush(mb.Compound):
    """ """
    def __init__(self, chain_length=4, alpha=pi/4):
        super(Brush, self).__init__()

        # Add parts
        self.add(Silane(), label='silane')
        self.add(Initiator(), label='initiator')
        self.add(mb.Polymer(MPC(alpha=alpha), n=chain_length,
                            port_labels=('up', 'down')), label='pmpc')
        self.add(CH3(), label='methyl')

        mb.force_overlap(self['initiator'], self['initiator']['down'], self['silane']['up'])
        mb.force_overlap(self['pmpc'], self['pmpc']['down'], self['initiator']['up'])
        mb.force_overlap(self['methyl'], self['methyl']['up'], self['pmpc']['up'])

        # Make self.port point to silane.bottom_port
        self.add(self['silane']['down'], label='down', containment=False)

if __name__ == "__main__":
    #pmpc = Brush(chain_length=1, alpha=pi/4)
    pmpc = Brush()
    print(pmpc)
