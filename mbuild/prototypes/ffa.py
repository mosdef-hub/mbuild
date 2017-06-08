import mbuild as mb
import numpy as np

from mbuild.lib.moieties.ch3 import CH3
from mbuild.prototypes.alkyl_monomer import AlkylMonomer
from mbuild.prototypes.cooh import COOH

class FFA(mb.Compound):
    """Creates a saturated free fatty acid of n carbons based on user
    input"""
    def __init__(self, chain_length, hcap=None):
        super(FFA, self).__init__()
        
        if hcap:
            self.add(COOH(), label='head')
        else:
            self.add(COOH(ester=True), label='head')
        
        tail = mb.Polymer(AlkylMonomer(), (chain_length - 2))
        self.add(tail, label='tail')
        mb.x_axis_transform(self['tail'], new_origin=self['tail'][0],
                point_on_x_axis=self['tail'][6],
                point_on_xy_plane=self['tail'][3])
        self.add(CH3(), label='tailcap')
        
        self['head'].rotate(np.pi, [0,1,0])
        self['head'].translate([-self['tail'][3].pos[0],
            self['tail'][3].pos[1], 0])
        mb.force_overlap(move_this=self['tailcap'],
                from_positions=self['tailcap']['up'],
                to_positions=self['tail']['up'])
        
        self['head']['up'].spin(-np.pi/2, self['head'][0].pos)
        mb.force_overlap(move_this=self['head'],
                from_positions=self['head']['up'],
                to_positions=self['tail']['down'])

if __name__ == '__main__':
    ffa = FFA(24)
    #ffa.energy_minimization() for single C-O bond character
    ffa.save('ffa24.mol2', overwrite=True)
