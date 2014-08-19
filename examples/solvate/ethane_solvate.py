import os
import sys
from examples.ethane.ethane import Ethane
from mbuild.box import Box
from mbuild.mol2file import load_mol2
from mbuild.tools import solvate

__author__ = 'sallai'


if __name__ == "__main__":
    ethane = Ethane()

    # box = ethane.boundingbox()
    #
    # print box

    host_box = Box(mins=[-9,-9,-9], maxes=[19,19,19])
    print("Host (ethane) box: {}".format(host_box))


    # Look for data file in same directory as this python module.
    current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
    new_path = os.path.join(current_dir, 'spc216.mol2')
    water = load_mol2(new_path)

    guest_box = water.boundingbox()

    print("Guest (water) box: {}".format(guest_box))

    # add water to ethane
    solvate(ethane, water, host_box, guest_box)

    from mbuild.plot import Plot
    Plot(ethane, verbose=False, atoms=True, bonds=True, angles=False, dihedrals=False).show()
