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

    host_box = Box(mins=[-9,-9,-9], maxes=[9,9,9])

    water_box = load_mol2("../../mbuild/components/spc216.mol2")

    guest_box = water_box.boundingbox()

    print guest_box


    solvate(ethane, water_box, host_box, guest_box)
    #
    #
    # from mbuild.plot import Plot
    # Plot(ethane, verbose=False, atoms=True, bonds=True, angles=False, dihedrals=False).show()
