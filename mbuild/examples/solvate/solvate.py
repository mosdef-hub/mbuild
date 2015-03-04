from __future__ import print_function

import mbuild as mb
from mbuild.testing.tools import get_fn
from mbuild.examples.ethane.ethane import Ethane


def main():
    """Solvate an ethane molecule in a Box of water. """

    # Create ethane and give it a box.
    ethane = Ethane()
    host_box = mb.Box(mins=[-.9, -.9, -.9], maxs=[1.9, 1.9, 1.9])
    print("Host (ethane) box: {}".format(host_box))

    # Create a water box.
    water = mb.Compound()
    water.append_from_file(get_fn('spc216.pdb'))
    guest_box = water.boundingbox()
    print("Guest (water) box: {}".format(guest_box))

    # Solvate ethane with water box.
    mb.solvate(ethane, water, host_box, guest_box)

    ethane = ethane.to_trajectory()
    ethane.topology.load_ff_bonds()
    ethane.save(filename='ethane.hoomdxml')

if __name__ == "__main__":
    main()
