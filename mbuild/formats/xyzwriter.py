import numpy as np

__all__ = ['write_xyz']


def write_xyz(structure, filename):
    """Output an XYZ file.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd structure object
    filename : str
        Path of the output file

    Notes
    -----
    Coordatates are written in Angstroms. This follows the convention for the
    XYZ file format.

    The XYZ file format neglects many important details, notably as bonds,
    residues, and box information.

    """

    xyz = np.array([[10.0*atom.xx, 10.0*atom.xy, 10.0*atom.xz] for atom in structure.atoms])
    types = [atom.name for atom in structure.atoms]

    with open(filename, 'w') as xyz_file:
        xyz_file.write(len(structure.atoms))
        xyz_file.write(filename+' - created by mBuild\n\n')
        for type, coords in zip(types, xyz):
            xyz_file.write('{:s}\t{:.3}\t{:.3f}\t{:.3f}\n'.format(type, *coords))
