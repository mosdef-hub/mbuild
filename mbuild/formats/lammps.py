from __future__ import print_function

import os
import operator

from mdtraj.utils import in_units_of
from mdtraj.formats.registry import _FormatRegistry
from six import string_types



__all__ = ['load_lammpsdata', 'LAMMPSTopologyFile']


@_FormatRegistry.register_loader('.lammps')
@_FormatRegistry.register_loader('.lmp')
def load_lammpsdata(filename, unitset=None):
    """Load a LAMMPS data file from disk.

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
    if not isinstance(filename, string_types):
        raise TypeError('Filename must be of type string for load_lammpstrj. '
                        'you supplied {0}'.format(type(filename)))

    with LAMMPSTopologyFile(filename, unitset=unitset) as f:
        compound = f.read()
    return compound


@_FormatRegistry.register_fileobject('.lammps')
@_FormatRegistry.register_fileobject('.lmp')
class LAMMPSTopologyFile(object):
    """ """
    per_particle_nodes = ['image', 'velocity', 'acceleration', 'mass',
                          'diameter', 'charge', 'body', 'orientation',
                          'moment_inertia']
    multi_particle_nodes = [('bond', 2), ('angle', 3), ('dihedral', 4),
                            ('improper', 4)]

    def __init__(self, filename, mode='r', force_overwrite=True, unitset='real'):
        """Open a LAMMPS data file for reading/writing. """

        self._is_open = False
        self._filename = filename
        self._mode = mode
        self._frame_index = 0
        # Track which line we're on. this is not essential, but its useful
        # when reporting errors to the user to say what line it occured on.
        # Looking at you LAMMPS...
        self._line_counter = 0

        if mode == 'r':
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            from mbuild.compound import Compound
            self._fh = open(filename, 'r')
            self._is_open = True
            self.compound = Compound()
        elif mode == 'w':
            if os.path.exists(filename) and not force_overwrite:
                raise IOError("The file '%s' already exists" % filename)
            self._fh = open(filename, 'w')
            self._is_open = True
        else:
            raise ValueError('mode must be one of "r" or "w". '
                             'you supplied "{0}"'.format(mode))

        self.u = dict()
        self.u['radians'] = 'radians'
        self.u['degrees'] = 'degrees'
        if unitset == 'real':
            self.u['distance'] = 'angstroms'
            self.u['velocity'] = 'angstroms/femtosecond'
            self.u['energy'] = 'kilocalories/mole'
            self.u['mass'] = 'grams/mole'
            self.u['charge'] = 'elementary_charge'
            self.u['mole'] = 'mole'
        else:
            raise Exception("Unsupported unit set specified: {0}".format(unitset))

    def read(self):
        """ """
        raise NotImplementedError

    def _read_unitcell_vectors(self, box):
        """Parse unitcell vectors from box node.  """
        pass

    def write(self, traj):
        """Output a Trajectory as a LAMMPS data file.

        Args:
            traj (md.Trajectory): The Trajectory to be output.
            filename (str, optional): Path of the output file.

        """
        directives_to_write = list()
        mass_list = list()
        mass_list.append('\n')
        mass_list.append('Masses\n')
        mass_list.append('\n')

        atom_list = list()
        atom_list.append('\n')
        atom_list.append('Atoms\n')
        atom_list.append('\n')

        numeric_types = dict()
        atom_type_n = 1
        for chain in traj.top.chains:
            for atom in chain.atoms:
                if (atom.name, atom.atomtype) not in numeric_types:
                    numeric_types[(atom.name, atom.atomtype)] = atom_type_n
                    mass = in_units_of(atom.element.mass, 'grams/moles', self.u['mass'])
                    mass_list.append('{0:d} {1:8.4f} # {2}, {3}\n'.format(
                        atom_type_n, mass, atom.name, atom.atomtype))
                    atom_type_n += 1
                x, y, z = in_units_of(traj.xyz[0][atom.index], 'nanometers', self.u['distance'])
                entry = '{0:-d} {1:-d} {2:-d} {3:5.8f} {4:8.5f} {5:8.5f} {6:8.5f}  # {7}, {8}\n'.format(
                    atom.index + 1, chain.index + 1, numeric_types[(atom.name, atom.atomtype)],
                    atom.charge, x, y, z, atom.name, atom.atomtype)
                atom_list.append(entry)

        n_atom_types = len(numeric_types)
        directives_to_write.append(mass_list)
        directives_to_write.append(atom_list)
        print('(string: numeric) types for atoms')
        print(sorted(numeric_types.items(), key=operator.itemgetter(1)))

        optional_directives = [('bonds', 2), ('angles', 3),
                               ('dihedrals', 4), ('impropers', 4)]
        number_of_terms = dict()
        number_of_types = dict()
        for directive, n_terms in optional_directives:
            number_of_terms[directive] = 0
            numeric_types = dict()
            if getattr(traj.top, '_ff_{0}'.format(directive)):
                list_of_terms = list()
                list_of_terms.append('\n')
                list_of_terms.append('{0}\n'.format(directive.title()))
                list_of_terms.append('\n')

                term_type_n = 1
                term_n = 0
                for term in getattr(traj.top, 'ff_{0}'.format(directive)):
                    if term.kind not in numeric_types:
                        numeric_types[term.kind] = term_type_n
                        term_type_n += 1
                    entry = '{0:-d} {1:d} '.format(term_n + 1, numeric_types[term.kind])
                    for n in range(n_terms):
                        entry += '{0:d} '.format(getattr(term, 'atom{0}'.format(n + 1)).index + 1)
                    entry += ' # {0}\n'.format(term.kind)
                    list_of_terms.append(entry)
                    term_n += 1
                number_of_terms[directive] = term_n
                directives_to_write.append(list_of_terms)
            number_of_types[directive] = len(numeric_types)
            print('(string: numeric types) for {0}'.format(directive))
            print(sorted(numeric_types.items(), key=operator.itemgetter(1)))

        self._fh.write('Generated by mBuild\n\n')

        n_bonds = number_of_terms['bonds']
        n_angles = number_of_terms['angles']
        n_dihedrals = number_of_terms['dihedrals']
        n_impropers = number_of_terms['impropers']

        n_bond_types = number_of_types['bonds']
        n_angle_types = number_of_types['angles']
        n_dihedral_types = number_of_types['dihedrals']
        n_improper_types = number_of_types['impropers']

        self._fh.write('{0} atoms\n'.format(traj.n_atoms))
        self._fh.write('{0} bonds\n'.format(n_bonds))
        self._fh.write('{0} angles\n'.format(n_angles))
        self._fh.write('{0} dihedrals\n'.format(n_dihedrals))
        self._fh.write('{0} impropers\n'.format(n_impropers))
        self._fh.write('\n')

        self._fh.write('{0} atom types\n'.format(n_atom_types))
        self._fh.write('{0} bond types\n'.format(n_bond_types))
        self._fh.write('{0} angle types\n'.format(n_angle_types))
        self._fh.write('{0} dihedral types\n'.format(n_dihedral_types))
        self._fh.write('{0} improper types\n'.format(n_improper_types))
        self._fh.write('\n')

        mins = traj.xyz[0].min(axis=0)
        maxs = mins + traj.unitcell_lengths[0]

        mins = in_units_of(mins, 'nanometers', self.u['distance'])
        maxs = in_units_of(maxs, 'nanometers', self.u['distance'])
        self._fh.write(
            '{0:10.6f} {1:10.6f} xlo xhi\n'.format(mins[0], maxs[0]))
        self._fh.write(
            '{0:10.6f} {1:10.6f} ylo yhi\n'.format(mins[1], maxs[1]))
        self._fh.write(
            '{0:10.6f} {1:10.6f} zlo zhi\n'.format(mins[2], maxs[2]))

        for directive in directives_to_write:
            for entry in directive:
                self._fh.write(entry)

    def close(self):
        """Close the LAMMPS data file. """
        if self._is_open:
            self._fh.close()
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
    import numpy as np

    from mbuild.examples.ethane.ethane import Ethane

    ethane = Ethane()
    ethane.save("ethane.lammps", forcefield='OPLS-AA')


