from __future__ import division, print_function

from numpy import pi

import mbuild as mb

from mbuild.examples.pmpc.brush import Brush


def main():
    """A quick test for updating a Compound and traversing its hierarchy. """
    # Create a compound and write it to file.
    brush1 = Brush()
    brush1.save("brush1.pdb")

    # Create another compound, rotate it and write it to file.
    brush2 = Brush()
    mb.rotate_around_z(brush2, pi/2)
    brush2.save("brush2.pdb")

    # Load brush2.pdb into brush1, modifying the atom positions of brush1.
    brush1.update_coordinates("brush2.pdb")
    brush1.save("modified_brush1.pdb")

    # Access the internals of the updated brush1.
    print(brush1.pmpc)

    for mpc in brush1.pmpc.monomer:
        print(mpc)

if __name__ == "__main__":
    main()
