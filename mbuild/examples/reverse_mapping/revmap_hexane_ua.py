import mbuild as mb
from mapping_moieties.ch2_aa import CH2_aa
from mapping_moieties.ch3_aa import CH3_aa

coarse_grained = mb.load('two_hexane_ua.mol2')
mapping_moieties = {'CH3':CH3_aa, 
                    'CH2':CH2_aa}

recovered = mb.reverse_map(coarse_grained, mapping_moieties)
