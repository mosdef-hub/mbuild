import mbuild as mb
import numpy as np
from mbuild.lib.UA_molecules.ch3ua import CH3UA
from mbuild.lib.UA_molecules.alkyl_monomerua import AlkylMonomerUA
from mbuild.lib.prototypes.OH import OH


class ALCUA(mb.Compound):
    """Creates a straight-chain alcohol of n carbons based on user
    input"""
    def __init__(self, chain_length):
        super(ALCUA, self).__init__()
        self.add(CH3UA(), label='tailcap')
        tail = mb.Polymer(AlkylMonomerUA(), (chain_length - 1))
        self.add(tail, label='tail')
        self.add(OH(), label='head')
        mb.force_overlap(move_this=self['tailcap'], from_positions=self['tailcap']['up'],
                         to_positions=self['tail']['down'])
        mb.force_overlap(move_this=self['head'], from_positions=self['head']['up'],
                         to_positions=self['tail']['up'])
        mb.z_axis_transform(self, new_origin=self['head'][0], point_on_z_axis=self['tailcap'][0])
        self.rotate(np.pi, [1, 0, 0])
        self.name = 'ALC' + str(chain_length)

if __name__ == '__main__':
    alcohol = ALCUA(16)
    alcohol.save('alcua.mol2', overwrite=True)
