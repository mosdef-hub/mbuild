import mbuild as mb
from mapping_moieties.ch2_aa import CH2_aa
from mapping_moieties.ch3_aa import CH3_aa

""" The coarse-grained system is two propanes,
each expressed with their united atom, implicit hydrogen
beads. The reverse mapping finds the CH3 and CH2 beads 
and replaces them with their all-atom representations"""

coarse_grained = mb.load('two_propane_ua.mol2')
mapping_moieties = {'_CH3':CH3_aa, 
                    '_CH2':CH2_aa}

recovered = mb.reverse_map(coarse_grained, mapping_moieties)
