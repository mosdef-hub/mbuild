import mbuild as mb
from mapping_moieties.propane_aa import Propane_aa

"""Here, the coarse-grained system is two hexanes, 
each expressed as two propanes. The reverse-mapping
finds the propanes and replaces them with
all-atom representations

"""
coarse_grained = mb.load('two_hexane_cg.mol2')
mapping_moieties = {'Propane': Propane_aa}

recovered = mb.reverse_map(coarse_grained, mapping_moieties, energy_minimize=True)
recovered.save('revmap.mol2',overwrite=True)
