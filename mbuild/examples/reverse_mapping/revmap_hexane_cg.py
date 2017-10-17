import mbuild as mb
from mapping_moieties.propane_aa import Propane_aa
from reverse_mapping import *

coarse_grained = mb.load('testing/hexane_cg.mol2')
mapping_moieties = {'Propane': Propane_aa}

recovered = mb.reverse_map(coarse_grained, mapping_moieties)

