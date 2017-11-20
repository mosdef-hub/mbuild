import mbuild as mb
from mapping_moieties.propane_ua import Propane_ua

"""Here, the coarse-grained system is two hexanes, 
each expressed as two propanes. The reverse-mapping
finds the propanes and replaces them with
united-atom representations

"""
coarse_grained = mb.load('two_hexane_cg.mol2')
mapping_moieties = {'Propane': Propane_ua}

recovered = mb.reverse_map(coarse_grained, mapping_moieties, minimize_energy=True,forcefield='toy_ua.xml')
recovered.save('revmap.mol2', overwrite=True)

