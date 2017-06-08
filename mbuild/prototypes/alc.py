import mbuild as mb
from mbuild.lib.moieties.ch3 import CH3
from mbuild.prototypes.alkyl_monomer import AlkylMonomer
from mbuild.prototypes.OH import OH

class ALC(mb.Compound):
    """Creates a straight-chain alcohol of n carbons based on user
    input"""
    def __init__(self, chain_length):
        super(ALC, self).__init__()
        self.add(CH3(), label='tailcap')
        tail = mb.Polymer(AlkylMonomer(), (chain_length - 1))
        self.add(tail, label='tail')
        self.add(OH(), label='head')
        mb.force_overlap(move_this=self['tailcap'],
                from_positions=self['tailcap']['up'],
                to_positions=self['tail']['up'])
        mb.force_overlap(move_this=self['head'],
                from_positions=self['head']['up'],
                to_positions=self['tail']['down'])

if __name__ == '__main__':
    alcohol = ALC(16)
    alcohol.energy_minimization()
    alcohol.save('OH16.mol2', overwrite=True)
