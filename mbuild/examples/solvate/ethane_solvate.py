from mbuild.examples.bilayer.eceramidens import ECeramideNS
from mbuild.formats.hoomdxml import save_hoomdxml

__author__ = 'sallai'

from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.testing.tools import get_fn
from mbuild.tools import solvate
from mbuild.examples.ethane.ethane import Ethane

if __name__ == "__main__":
    ethane = Ethane()

    host_box = Box(mins=[-.9, -.9, -.9], maxs=[1.9, 1.9, 1.9])
    print("Host (ethane) box: {}".format(host_box))

    water = Compound()
    water.append_from_file(get_fn('spc216.pdb'))
    #water.append_from_file(get_fn('water.hoomdxml'))

    guest_box = water.boundingbox()
    print("Guest (water) box: {}".format(guest_box))

    # add water to ethane
    solvate(ethane, water, host_box, guest_box)

    ethane = ethane.to_trajectory()
    ethane.topology.load_ff_bonds()
    save_hoomdxml(ethane, filename='ethane.hoomdxml')
    #from mbuild.plot import Plot
    #Plot(ethane, verbose=False, atoms=True, bonds=True, angles=False, dihedrals=False).show()
