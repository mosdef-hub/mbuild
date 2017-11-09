import mbuild as mb

from mbuild.lib.prototypes.pctails import PCTails
from mbuild.lib.prototypes.pchead import PCHead


class DMPC(mb.Compound):
    def __init__(self):
        super(DMPC, self).__init__()
        
        self.add(PCHead(), label='headgroup')
        self.add(PCTails(14, 14), label='ffatails')
        
        mb.force_overlap(move_this=self['ffatails'], from_positions=self['ffatails']['base']['up'],
                         to_positions=self['headgroup']['alkyl_split']['up'])

if __name__ == '__main__':
    dmpc = DMPC()
    dmpc.save('dmpc.mol2', overwrite=True)
