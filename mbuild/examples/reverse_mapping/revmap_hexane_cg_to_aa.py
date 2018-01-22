import mbuild as mb
from mbuild.examples import Alkane

"""Here, the coarse-grained system is two hexanes, 
each expressed as two propanes. The reverse-mapping
finds the propanes and replaces them with
all-atom representations

"""
coarse_grained = mb.load('one_hexane_cg.mol2')
mapping_moieties = {'Propane': Alkane(n=3, cap_end=False)}

recovered = mb.reverse_map(coarse_grained, mapping_moieties, energy_minimize=True,forcefield='oplsaa.xml' ,scale_torsions=0.10)
recovered.save('revmap.mol2',overwrite=True)
