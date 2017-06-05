import mbuild as mb
from mbuild.lib.moieties.ch3 import CH3
from mbuild.prototypes.alkyl_monomer import AlkylMonomer
from mbuild.prototypes.cooh import COOH

class FFA(mb.Compound):
    """Creates a saturated free fatty acid of n carbons based on user
    input"""
    def __init__(self, chain_length):
        super(FFA, self).__init__()
        self.add(CH3(), label='tailcap')
        tail = mb.Polymer(AlkylMonomer(), (chain_length - 2))
        self.add(tail, label='tail')
        self.add(COOH(), label='head')
        mb.force_overlap(move_this=self['tailcap'],
                from_positions=self['tailcap']['up'],
                to_positions=self['tail']['up'])
        mb.force_overlap(move_this=self['head'],
                from_positions=self['head']['down'],
                to_positions=self['tail']['down'])

if __name__ == '__main__':
    ffa = FFA(16)
    ffa.energy_minimization()
    ffa.save('ffa16.mol2', overwrite=True)
