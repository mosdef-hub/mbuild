import numpy as np

import mbuild as mb

__all__ = ['read_xyz']


def read_xyz(filename):
    """Read an XYZ file.

    Parameters
    ----------
    filename : str
        Path of the output file

    Returns
    -------
    compound : mb.Compound

    Notes
    -----
    The XYZ file format neglects many important details, notably as bonds,
    residues, and box information.

    """

    compound = mb.Compound()

    with open(filename, 'r') as xyz_file:
        n_atoms = int(xyz_file.readline())
        xyz_file.readline()
        coords = np.zeros(shape=(n_atoms, 3), dtype=np.float64)
        for row, _ in enumerate(coords):
            line = xyz_file.readline().split()
            coords[row] = line[1:4]
            coords[row] *= 0.1
            particle = mb.Compound(pos=coords[row], name=line[0])
            compound.add(mb.clone(particle))

    return compound
