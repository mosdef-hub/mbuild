import mbuild as mb
import numpy as np
from mbuild.lib.moieties.ch3 import CH3
from mbuild.lib.prototypes.alkyl_monomer import AlkylMonomer
from mbuild.lib.prototypes.OH import OH


class ALC(mb.Compound):
    """Creates a straight-chain alcohol of n carbons based on user
    input"""
    def __init__(self, chain_length):
        super(ALC, self).__init__()
        self.add(CH3(), label='tailcap')
        tail = mb.Polymer(AlkylMonomer(), (chain_length - 1))
        self.add(tail, label='tail')
        self.add(OH(), label='head')
        mb.force_overlap(move_this=self['tailcap'], from_positions=self['tailcap']['up'],
                         to_positions=self['tail']['up'])
        mb.force_overlap(move_this=self['head'], from_positions=self['head']['up'],
                         to_positions=self['tail']['down'])
        mb.z_axis_transform(self, new_origin=self['head'][0], point_on_z_axis=self['tailcap'][0])
        self.rotate(np.pi, [1, 0, 0])

if __name__ == '__main__':
    alcohol = ALC(16)
    alcohol.save('alc.mol2', overwrite=True)
