import mbuild as mb
from mapping_moieties.ch2_aa import CH2_aa
from mapping_moieties.ch3_aa import CH3_aa

""" The coarse grained system is two hexanes,
each expressed as their united atom models with 
implicit hydrogens. The reverse-mapping finds 
all the united atom beads (CH3 and CH2) and replaces
them with their all-atom representtations"""
coarse_grained = mb.load('new.mol2')
#coarse_grained = mb.load('two_hexane_ua.mol2')
mapping_moieties = {'_CH3':CH3_aa(), 
                    '_CH2':CH2_aa()}

recovered = mb.reverse_map(coarse_grained, mapping_moieties)
recovered.save('revmap.mol2',overwrite=True)
