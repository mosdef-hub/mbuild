from hoomd_script import *

init.read_xml(filename='ethane.xml')

harmonic = bond.harmonic()
harmonic.set_coeff('C-C', k=300, r0=0.8)
harmonic.set_coeff('C-H', k=150, r0=0.6)

lj = pair.lj(r_cut=3.0)
lj.pair_coeff.set('C', 'C', epsilon=1.0, sigma=2.0)
lj.pair_coeff.set('H', 'C', epsilon=1.0, sigma=1.0, alpha=0.5)
lj.pair_coeff.set('H', 'H', epsilon=1.0, sigma=1.0)

ethane = group.ethane()
fire = integrate.mod_minimize_fire(group=group.all(), dt=0.05, ftol=1e-2, etol=1e-7)

integrate.mode_standard(dt=0.005)
integrate.nvt(group=ethane, T=1, tau=0.5)

dump.dcd(filename='ethane.dcd', period=100)
run(1e3)
