import mbuild as mb
import numpy as np

from mbuild.lib.prototypes.OH import OH


class COOH(mb.Compound):
    """Creates headgroup of a carboxylic acid"""
    def __init__(self, ester=None, cis=False):
        super(COOH, self).__init__()
        if ester:
            self.add(mb.Particle(name='OE'), label='O[$]')
            self.add(mb.Particle(name='C'), label='C')
            self.add(mb.Particle(name='O', pos=[0, .123, 0]), label='O[$]')
            self.add_bond((self['C'], self['O'][1]))

            self.add(mb.Port(anchor=self['C'], orientation=[-1, -np.tan(50.5*np.pi/180), 0],
                             separation=.132/2), label='up')

            self['O'][0].translate([.15, 0, 0])
            theta = ((180-111) / 2) * np.pi / 180
            self['O'][0].rotate(-theta, around=[0, 0, 1])

            self.add_bond((self['C'], self['O'][0]))
            
            if cis:
                self.add(mb.Port(anchor=self['O'][0], orientation=[0, -1, 0], separation=.14/2),
                         label='down')
            else:
                self.add(mb.Port(anchor=self['O'][0], orientation=[1, np.tan(52.5*np.pi/180), 0],
                                 separation=.14/2), label='down')
        else:
            self.add(mb.Particle(name='O', pos=[0, .123, 0]), label='O[$]')
            self.add(mb.Particle(name='C'), label='C')
            self.add(OH(), label='OH')
            self.add_bond((self['C'], self['O'][0]))
            self.add(mb.Port(anchor=self['C'], orientation=[1, -np.tan(32*np.pi/180), 0],
                             separation=.132/2), label='down')

            self.add(mb.Port(anchor=self['C'], orientation=[-1, -np.tan(34.5*np.pi/180), 0],
                             separation=.132/2), label='up')
            
            mb.force_overlap(move_this=self['OH'], from_positions=self['OH']['up'], to_positions=self['down'])

if __name__ == '__main__':
    cooh = COOH(ester=True)
    cooh.save('cooh.mol2', overwrite=True)
