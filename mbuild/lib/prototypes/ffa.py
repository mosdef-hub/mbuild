import mbuild as mb
import numpy as np

from mbuild.lib.prototypes import COOH, AlkylMonomer
from mbuild.lib.moieties.ch3 import CH3


class FFA(mb.Compound):
    """Creates a saturated free fatty acid of n carbons based on user
    input for the united atom model"""
    def __init__(self, chain_length, ester=True):
        super(FFA, self).__init__()
        
        if ester:
            self.add(COOH(ester=True), label='head')
        
        tail = mb.Polymer(AlkylMonomer(), (chain_length - 2))
        if not ester:
            self.add(CH3(), label='tailcap')
        self.add(tail, label='tail')
        if ester:
            self.add(CH3(), label='tailcap')
        mb.x_axis_transform(self['tail'], new_origin=self['tail'][0], point_on_x_axis=self['tail'][2],
                            point_on_xy_plane=self['tail'][1])

        if not ester:
            self.add(COOH(), label='head')
        if ester:
            mb.force_overlap(move_this=self['tailcap'], from_positions=self['tailcap']['up'],
                             to_positions=self['tail']['up'])
        else:
            mb.force_overlap(move_this=self['tailcap'], from_positions=self['tailcap']['up'],
                             to_positions=self['tail']['down'])

        self['head']['up'].spin(-np.pi/2, self['head']['up'].pos)
        if ester: 
            mb.force_overlap(move_this=self['head'], from_positions=self['head']['up'],
                             to_positions=self['tail']['down'])
            self.spin(np.pi/2, [0, 1, 0])
        else:
            mb.force_overlap(move_this=self['head'], from_positions=self['head']['up'],
                             to_positions=self['tail']['up'])
            self.spin(-np.pi/2, [0, 1, 0])
        self.name = 'FFA' + str(chain_length)

if __name__ == '__main__':
    ffa = FFA(24)
    ffa.save('ffa.mol2', overwrite=True)
