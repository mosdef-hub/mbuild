from __future__ import print_function

import os

#from mdtraj.utils import in_units_of
from mdtraj.formats.registry import _FormatRegistry
from six import string_types

__all__ = ['load_gromacs', 'GROMACSTopologyFile']


@_FormatRegistry.register_loader('.gro')
@_FormatRegistry.register_loader('.top')
def load_gromacstop(top_filename=None, gro_filename=None):
    """Load a GROMACS .top file and accompanying .gro from disk.

    Parameters
    ----------
    filename : str
        Path to data file.
    optional_nodes : str, optional, default=['bonds']
        Read specified nodes in file other than 'box', 'position' and 'type'.

    Returns
    -------
    compound : mb.Compound

    """
    # TODO: If only one is given, assume the other has the same name, different
    #  extension.

    # if not isinstance(filename, string_types):
    #     raise TypeError('Filename must be of type string for load_lammpstrj. '
    #                     'you supplied {0}'.format(type(filename)))
    #
    # with GROMACSTopologyFile(filename) as f:
    #     compound = f.read()
    # return compound
    raise NotImplementedError


@_FormatRegistry.register_fileobject('.gro')
@_FormatRegistry.register_fileobject('.top')
class GROMACSTopologyFile(object):
    """ """
    def __init__(self, top_filename, gro_filename, mode='r', force_overwrite=True):
        """Open a GROMACS .gro/.top file pair for reading/writing. """

        self._is_open = False
        self._top_filename = top_filename
        self._gro_filename = gro_filename
        self._mode = mode
        self._frame_index = 0
        # Track which line we're on. this is not essential, but its useful
        # when reporting errors to the user to say what line it occured on.
        # Looking at you GROMACS...
        self._line_counter = 0

        if mode == 'r':
            if not os.path.exists(top_filename):
                raise IOError("The file '{0}' doesn't exist".format(top_filename))
            if not os.path.exists(gro_filename):
                raise IOError("The file '{0}' doesn't exist".format(gro_filename))
            self._topfh = open(top_filename, 'r')
            self._grofh = open(gro_filename, 'r')
            self._is_open = True

            from mbuild.compound import Compound
            self.compound = Compound()
        elif mode == 'w':
            if os.path.exists(top_filename) and not force_overwrite:
                raise IOError("The file '{0}' already exists".format(top_filename))
            if os.path.exists(gro_filename) and not force_overwrite:
                raise IOError("The file '{0}' already exists".format(gro_filename))
            self._topfh = open(top_filename, 'w')
            self._grofh = open(gro_filename, 'w')
            self._is_open = True
        else:
            raise ValueError('mode must be one of "r" or "w". '
                             'you supplied "{0}"'.format(mode))

        self.u = dict()
        self.u['radians'] = 'radians'
        self.u['degrees'] = 'degrees'
        self.u['distance'] = 'nanometers'
        self.u['velocity'] = 'nanometers/picosecond'
        self.u['energy'] = 'kilojoules/mole'
        self.u['mass'] = 'grams/mole'
        self.u['charge'] = 'elementary_charge'
        self.u['mole'] = 'mole'

    def read(self):
        """ """
        raise NotImplementedError

    def _read_unitcell_vectors(self, box):
        """Parse unitcell vectors from box node.  """
        pass

    def write(self, traj):
        """Output a Trajectory as a GROMACS data file.

        Args:
            traj (md.Trajectory): The Trajectory to be output.
            filename (str, optional): Path of the output file.

        """
        basename = os.path.splitext(self._top_filename)[0]
        # .top file
        if traj.top.forcefield == 'opls-aa':
            self._topfh.write('#include "oplsaa.ff/forcefield.itp"\n\n')

        # TODO: iteration over chains
        self._topfh.write('\n[ moleculetype ]\n')
        self._topfh.write('; name nrexcl\n')
        self._topfh.write('mbuild 3\n')

        self._topfh.write('\n[ atoms ]\n')
        for i, atom in enumerate(traj.topology._atoms):
            self._topfh.write('{:8d} {:5s} {:8d} {:5s} {:5s} {:8d} {:8.4f} {:8.4f}\n'.format(
                atom.index + 1, atom.atomtype, atom.residue.index,
                atom.residue.name,
                atom.bondtype, i+1, atom.charge, atom.element.mass))

        if traj.topology._ff_pairs:
            self._topfh.write('\n[ pairs ]\n')
            for atom1, atom2 in traj.topology.ff_pairs:
                self._topfh.write('{:d} {:d}\n'.format(
                    atom1.index + 1, atom2.index + 1))

        if traj.topology._ff_bonds:
            self._topfh.write('\n[ bonds ]\n')
            for bond in traj.topology.ff_bonds:
                self._topfh.write('{:d} {:d} {:d}\n'.format(
                    bond.atom1.index + 1, bond.atom2.index + 1, 1))

        if traj.topology._ff_angles:
            self._topfh.write('\n[ angles ]\n')
            for angle in traj.topology.ff_angles:
                self._topfh.write('{:d} {:d} {:d} {:d}\n'.format(
                    angle.atom1.index + 1, angle.atom2.index + 1,
                    angle.atom3.index + 1, 1))

        if traj.topology._ff_dihedrals:
            self._topfh.write('\n[ dihedrals ]\n')
            self._topfh.write('; impropers\n')
            for dihedral in traj.topology.ff_dihedrals:
                self._topfh.write('{:d} {:d} {:d} {:d} {:d}\n'.format(
                    dihedral.atom1.index + 1, dihedral.atom2.index + 1,
                    dihedral.atom3.index + 1, dihedral.atom4.index + 1, 3))

        self._topfh.write('\n[ system ]\n')
        self._topfh.write('{}\n'.format(basename))

        self._topfh.write('\n[ molecules ]\n')
        self._topfh.write('mbuild 1\n')

        # .gro file
        self._grofh.write('{}\n'.format(basename))
        self._grofh.write('{}\n'.format(traj.n_atoms))
        for atom, xyz in zip(traj.top._atoms, traj.xyz[0]):
            self._grofh.write('{:5d}{:<5s}{:5s}{:5d}'
                    '{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}\n'.format(
                atom.residue.index, atom.residue.name, atom.bondtype, atom.index + 1,
                xyz[0], xyz[1], xyz[2], 0.0, 0.0, 0.0))
        box = traj.unitcell_vectors[0]
        self._grofh.write('{:10.5f}{:10.5f}{:10.5f}'
                          '{:10.5f}{:10.5f}{:10.5f}'
                          '{:10.5f}{:10.5f}{:10.5f}'.format(
            box[0, 0], box[1, 1], box[2, 2],
            box[1, 0], box[2, 0], box[0, 1],
            box[2, 1], box[0, 2], box[1, 2]))
        self._grofh.write('\n')

    def close(self):
        """Close the GROMACS file. """
        if self._is_open:
            self._topfh.close()
            self._grofh.close()
            self._is_open = False

    def __del__(self):
        self.close()

    def __enter__(self):
        """Support the context manager protocol. """
        return self

    def __exit__(self, *exc_info):
        """Support the context manager protocol. """
        self.close()


if __name__ == "__main__":
    from mbuild.examples.ethane.ethane import Ethane

    Ethane().save('ethane.gro', forcefield='opls-aa')

