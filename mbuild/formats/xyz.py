import numpy as np

import mbuild as mb
from mbuild.exceptions import MBuildError
from ele import element_from_symbol
from ele.exceptions import ElementError
from warnings import warn

__all__ = ['read_xyz', 'write_xyz']


def read_xyz(filename, compound=None):
    """Read an XYZ file. The expected format is as follows:
    The first line contains the number of atoms in the file The second line
    contains a comment, which is not read.  Remaining lines, one for each
    atom in the file, include an elemental symbol followed by X, Y, and Z
    coordinates in Angstroms. Columns are expected tbe separated by
    whitespace. See https://openbabel.org/wiki/XYZ_(format).

    Parameters
    ----------
    filename : str
        Path of the input file

    Returns
    -------
    compound : mb.Compound

    Notes
    -----
    The XYZ file format neglects many important details, including bonds,
    residues, and box information.

    There are some other flavors of the XYZ file format and not all are
    guaranteed to be compatible with this reader. For example, the TINKER
    XYZ format is not expected to be properly read.
    """

    if compound is None:
        compound = mb.Compound()

    guessed_elements = set()

    with open(filename, 'r') as xyz_file:
        n_atoms = int(xyz_file.readline())
        xyz_file.readline()
        coords = np.zeros(shape=(n_atoms, 3), dtype=np.float64)
        for row, _ in enumerate(coords):
            line = xyz_file.readline().split()
            if not line:
                msg = ('Incorrect number of lines in input file. Based on the '
                       'number in the first line of the file, {} rows of atoms '
                       'were expected, but at least one fewer was found.')
                raise MBuildError(msg.format(n_atoms))
            coords[row] = line[1:4]
            coords[row] *= 0.1
            name = line[0]
            try:
                element = name.capitalize()
                element = element_from_symbol(element)
            except ElementError:
                if name not in guessed_elements:
                    warn(
                        "No matching element found for {}; the particle will "
                        "be added to the compound without an element "
                        "attribute.".format(name)
                    )
                    guessed_elements.add(name)
                element = None

            particle = mb.Compound(pos=coords[row], name=name, element=element)
            compound.add(particle)

        # Verify we have read the last line by ensuring the next line in blank
        line = xyz_file.readline().split()
        if line:
            msg = ('Incorrect number of lines in input file. Based on the '
                   'number in the first line of the file, {} rows of atoms '
                   'were expected, but at least one more was found.')
            raise MBuildError(msg.format(n_atoms))

    return compound


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

    if isinstance(structure, mb.Compound):
        raise ValueError(
            'Expected a ParmEd structure, got an mbuild.Compound'
        )

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])
    types = [atom.name for atom in structure.atoms]

    with open(filename, 'w') as xyz_file:
        xyz_file.write(str(len(structure.atoms)))
        xyz_file.write('\n' + filename+' - created by mBuild\n')
        for typ, coords in zip(types, xyz):
            xyz_file.write('{:s} {:11.6f} {:11.6f} {:11.6f}\n'.format(typ, *coords))
