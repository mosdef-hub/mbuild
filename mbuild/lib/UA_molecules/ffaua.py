import mbuild as mb
import numpy as np

from mbuild.lib.UA_molecules.ch3ua import CH3UA
from mbuild.lib.UA_molecules.alkyl_monomerua import AlkylMonomerUA
from mbuild.lib.prototypes.cooh import COOH


class FFAUA(mb.Compound):
    """Creates a saturated free fatty acid of n carbons based on user
    input for the united atom model"""
    def __init__(self, chain_length, ester=True):
        super(FFAUA, self).__init__()
        
        if ester:
            self.add(COOH(ester=True), label='head')
        
        tail = mb.Polymer(AlkylMonomerUA(), (chain_length - 2))
        if not ester:
            self.add(CH3UA(), label='tailcap')
        self.add(tail, label='tail')
        if ester:
            self.add(CH3UA(), label='tailcap')
        mb.x_axis_transform(self['tail'], new_origin=self['tail'][0], point_on_x_axis=self['tail'][2],
                            point_on_xy_plane=self['tail'][1])

        if not ester:
            self.add(COOH(), label='head')
        self['head'].rotate(np.pi, [0, 1, 0])
        self['head'].translate([-self['tail'][1].pos[0], self['tail'][1].pos[1], 0])
        if ester:
            mb.force_overlap(move_this=self['tailcap'], from_positions=self['tailcap']['up'],
                             to_positions=self['tail']['up'])
        else:
            mb.force_overlap(move_this=self['tailcap'], from_positions=self['tailcap']['up'],
                             to_positions=self['tail']['down'])

        self['head']['up'].spin(np.pi/2, self['head']['C'].pos)
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
    ffa = FFAUA(24, ester=False)
    ffa.save('ffaua.mol2', overwrite=True)
