from __future__ import print_function

import mbuild as mb
from mbuild.components.small_groups.h2o import H2O
from mbuild.examples.ethane.ethane import Ethane


def main():
    """Solvate an ethane molecule in a Box of water. """

    # Create ethane and give it a box.
    ethane = Ethane()
    host_box = mb.Box(mins=[-1, -1, -1], maxs=[1, 1, 1])
    print("Host box: {}".format(host_box))

    # Create a water box.
    water = H2O()

    # Solvate ethane with water box.
    return mb.solvate(ethane, water, 500, host_box)


if __name__ == "__main__":
    solvated_ethane = main()
    #solvated_ethane.visualize()

