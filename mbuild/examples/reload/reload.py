import numpy as np

from mbuild.examples.pmpc_brush_layer.brush import Brush
from mbuild.coordinate_transform import rotate_around_z


__author__ = 'sallai'

def main():
    # create a compound and write it to file
    brush1 = Brush()
    brush1.save("brush1.pdb")

    # create another compound, rotate it and write it to file
    brush2 = Brush()
    rotate_around_z(brush2, np.pi/2)
    brush2.save("brush2.pdb")

    # load brush2.pdb into brush1, modifying the atom positions of brush1
    brush1.update_from_file("brush2.pdb")
    brush1.save("modified_brush1.pdb")

    # access the internals of brush2
    print brush1.pmpc

    for mpc in brush1.pmpc.monomer:
        print mpc

if __name__ == "__main__":
    main()
