import mbuild as mb
import numpy as np

from mbuild.UA_molecules.ch3ua import CH3UA
from mbuild.UA_molecules.alkyl_monomerua import AlkylMonomerUA
from mbuild.prototypes.cooh import COOH

class FFAUA(mb.Compound):
    """Creates a saturated free fatty acid of n carbons based on user
    input for the united atom model"""
    def __init__(self, chain_length, hcap=None):
        super(FFAUA, self).__init__()
        
        if hcap:
            self.add(COOH(), label='head')
        else:
            self.add(COOH(ester=True), label='head')
        
        tail = mb.Polymer(AlkylMonomerUA(), (chain_length - 2))
        self.add(tail, label='tail')
        mb.x_axis_transform(self['tail'], new_origin=self['tail'][0],
                point_on_x_axis=self['tail'][6],
                point_on_xy_plane=self['tail'][3])
        self.add(CH3UA(), label='tailcap')
        
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
        self.spin(np.pi/2, [0,1,0])

if __name__ == '__main__':
    ffa = FFAUA(24, hcap=True)
    #ffa.energy_minimization() for single C-O bond character
    ffa.save('ffa24ua.mol2', overwrite=True)
