# -*- coding: utf-8 -*-


# -- ==cg_hexane== --
import mbuild as mb

from mbuild.examples.alkane.alkane import Alkane


class Propane(mb.Compound):
    def __init__(self):
        super(Propane, self).__init__()

        c = Alkane(n=3, cap_front=True, cap_end=False)
        self.add(c, 'propane')

        self.add(c['down'], 'down', containment=False)

class Hexane(mb.Compound):
    def __init__(self):
        super(Hexane, self).__init__()

        self.add(Propane(), 'propane1')
        self.add(Propane(), 'propane2')

        mb.force_overlap(self['propane1'], 
                         self['propane1']['down'], 
                         self['propane2']['down'])

# -- ==cg_hexane== --