import mbuild as mb
from mbuild.lib.moieties import CH3
from numpy import pi
import parmed as pmd
from parmed.periodic_table import AtomicNum

from foyer.forcefield import apply_forcefield

from mbuild.examples.pmpc.mpc import MPC

chain = mb.Polymer(MPC(alpha=pi/4), n=50, port_labels=('up', 'down'))

chain.add(CH3(), 'methyl_front')
mb.equivalence_transform(chain['methyl_front'], chain['methyl_front']['up'], chain['up'])

chain.add(CH3(), 'methyl_end')
mb.equivalence_transform(chain['methyl_end'], chain['methyl_end']['up'], chain['down'])

struc = chain.to_parmed(title='pmpc')

#struc = pmd.load_file('mpc.mol2', structure=True)

for atom in struc.atoms:
    atom.element = AtomicNum[atom.name[:1]]

pmpc = apply_forcefield(struc, 'oplsaa', debug=False)
pmpc.save('pmpc.top', overwrite=True)
pmpc.save('pmpc.gro', overwrite=True)




