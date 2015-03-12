import os
from warnings import warn

from mdtraj.formats.registry import _FormatRegistry
import numpy as np
from six import string_types
from xml.etree import cElementTree

import mbuild.compound


__all__ = ['load_hoomdxml', 'HOOMDTopologyFile']


@_FormatRegistry.register_loader('.hoomdxml')
def load_hoomdxml(filename, lj_units=None):
    """Load a HOOMD-blue XML file form disk.

    Note: lj_units need to be normalized by nm, kJ/mol, and amu

    Required nodes for valid HOOMD simulation: box, position and type.

    Parameters
    ----------
    filename : str
        Path to xml file.

    Returns
    -------
    compound : mb.Compound

    """
    if not isinstance(filename, string_types):
        raise TypeError('Filename must be of type string for load_lammpstrj. '
                        'you supplied {0}'.format(type(filename)))

    with HOOMDTopologyFile(filename, lj_units=lj_units) as f:
        compound = f.read()
    return compound


@_FormatRegistry.register_fileobject('.hoomdxml')
class HOOMDTopologyFile(object):
    """ """
    per_particle_nodes = ['image', 'velocity', 'acceleration', 'mass',
                          'diameter', 'charge', 'body', 'orientation',
                          'moment_inertia']
    multi_particle_nodes = [('bond', 2), ('angle', 3), ('dihedral', 4),
                            ('improper', 4)]

    def __init__(self, filename, mode='r', force_overwrite=True, lj_units=None):
        """Open a HOOMD xml file for reading/writing. """

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
            self._fh = open(filename, 'r')
            self._is_open = True
            self.compound = mbuild.compound.Compound()
        elif mode == 'w':
            if os.path.exists(filename) and not force_overwrite:
                raise IOError("The file '%s' already exists" % filename)
            self._fh = open(filename, 'w')
            self._is_open = True
        else:
            raise ValueError('mode must be one of "r" or "w". '
                             'you supplied "{0}"'.format(mode))

        # Fundamental LJ units.
        if lj_units is None:
            self.u = {'distance': 1.0,
                      'energy': 1.0,
                      'mass': 1.0}
        else:
            assert isinstance(lj_units, dict)
            assert 'distance' in lj_units
            assert 'energy' in lj_units
            assert 'mass' in lj_units
            self.u = lj_units

        # Other derived LJ units.
        self.u['time'] = (np.sqrt(self.u['mass'] * self.u['distance']**2.0 / self.u['energy']))
        self.u['velocity'] = self.u['distance'] / self.u['time']
        self.u['acceleration'] = self.u['distance'] / self.u['time']**2.0
        self.u['diameter'] = self.u['distance']
        self.u['charge'] = 1.0
        # TODO: figure out charge
        self.u['moment_inertia'] = self.u['mass'] * self.u['distance']**2.0
        self.u['image'] = 1.0
        self.u['body'] = 1.0
        self.u['orientation'] = 1.0

    def read(self):
        """ """

        tree = cElementTree.parse(self._filename)
        self._config = tree.getroot().find('configuration')

        # Unitcell info.
        box = self._config.find('box')
        unitcell_vectors = self._read_unitcell_vectors(box) * self.u['distance']
        self.compound.periodicity = np.diag(unitcell_vectors)

        # Coordinates and atom types.
        coords = self._config.find('position').text.splitlines()[1:]
        types = self._config.find('type').text.splitlines()[1:]
        n_atoms = len(coords)
        atom_mapping = dict()
        for n, type_coord in enumerate(zip(types, coords)):
            atom_type, xyz = type_coord
            new_atom = mbuild.compound.Atom(str(atom_type), [float(x) * self.u['distance']
                                             for x in xyz.split()])
            self.compound.add(new_atom, label="{0}[$]".format(new_atom.kind))
            atom_mapping[n] = new_atom

        self._read_per_particle_nodes()
        self._read_multi_particle_nodes(atom_mapping)

        # TODO: read wall
        # <wall> has its information stored as attributes
        return self.compound

    def _read_per_particle_nodes(self):
        """ """
        for node in self.per_particle_nodes:
            try:
                node_text = self._config.find(node).text.splitlines()[1:]
            except AttributeError:
                warn("Specified node '{0}' does not exist in {1}".format(
                    node, self._filename))
            else:
                for raw_line, atom in zip(node_text, self.compound.atoms):
                    # TODO: not robust when e.g. charges are provided as ints
                    parsed_line = [int(x) if x.isdigit() else float(x) * self.u[node]
                                   for x in raw_line.split()]
                    if len(parsed_line) == 1:
                        parsed_line = parsed_line[0]
                    atom.extras[node] = parsed_line

    def _read_multi_particle_nodes(self, atom_mapping):
        """ """

        for node, n_indices in self.multi_particle_nodes:
            try:
                node_text = self._config.find(node).text.splitlines()[1:]
            except AttributeError:
                # Node doesn't exist.
                pass
            else:
                if node != 'bond':
                    self.compound.extras[node] = list()
                for raw_line in node_text:
                    parsed_line = [int(x) if x.isdigit() else x for x in raw_line.split()]
                    if node == 'bond':
                        atom1 = atom_mapping[parsed_line[1]]
                        atom2 = atom_mapping[parsed_line[2]]
                        new_bond = mbuild.compound.Bond(atom1, atom2)
                        self.compound.add(new_bond)
                    else:
                        self.compound.extras[node] = list()

    def _read_unitcell_vectors(self, box):
        """Parse unitcell vectors from box node.  """

        # TODO: make less horrible.
        for L in ['lx', 'LX', 'lX', 'Lx']:
            try:
                lx = float(box.attrib[L])
                break
            except KeyError:
                pass
        else:
            raise ValueError('Unable to find box length in x direction')

        for L in ['ly', 'LY', 'lY', 'Ly']:
            try:
                ly = float(box.attrib[L])
                break
            except KeyError:
                pass
        else:
            raise ValueError('Unable to find box length in y direction')

        for L in ['lz', 'LZ', 'lZ', 'Lz']:
            try:
                lz = float(box.attrib[L])
                break
            except KeyError:
                pass
        else:
            raise ValueError('Unable to find box length in z direction')

        try:
            xy = float(box.attrib['xy'])
            xz = float(box.attrib['xz'])
            yz = float(box.attrib['yz'])
        except KeyError:
            xy = 0.0
            xz = 0.0
            yz = 0.0

        return np.array([[lx,   xy*ly, xz*lz],
                         [0.0,     ly, yz*lz],
                         [0.0,   0.0,    lz]])

    def write(self, traj):
        """Output a Trajectory as a HOOMD XML file.

        Args:
            traj (md.Trajectory): The Trajectory to be output.
            filename (str, optional): Path of the output file.

        """
        self._fh.write("""<?xml version="1.3" encoding="UTF-8"?>\n""")
        self._fh.write("""<hoomd_xml>\n""")
        self._fh.write("""<configuration time_step="0">\n""")

        lx, ly, lz = traj.unitcell_lengths[0]
        xy = 0
        xz = 0
        yz = 0

        self._fh.write("""<box lx="{0}" ly="{1}" lz="{2}" xy="{3}" xz="{4}" yz="{5}" />\n""".format(
            lx, ly, lz, xy, xz, yz))

        self._fh.write("""<position num="{0}">\n""".format(traj.n_atoms))
        for x, y, z in traj.xyz[0]:
            self._fh.write("{0:8.5f} {1:8.5f} {2:8.5f}\n".format(x, y, z))
        self._fh.write("</position>\n")

        self._fh.write("""<type num="{0}">\n""".format(traj.n_atoms))
        for atom in traj.top.atoms:
            try:
                atomtype = atom.atomtype
            except AttributeError:
                atomtype = atom.name
            self._fh.write("{0}\n".format(atomtype))
        self._fh.write("</type>\n")

        optional_directives = [('bonds', 2), ('angles', 3), ('dihedrals', 4),
                               ('impropers', 4)]
        for directive, n_terms in optional_directives:
            if getattr(traj.top, '_ff_{0}'.format(directive)):
                self._fh.write("""<{0}>\n""".format(directive[:-1]))
                for term in getattr(traj.top, 'ff_{0}'.format(directive)):
                    entry = '{0}'.format(term.kind)
                    for n in range(n_terms):
                        entry += ' {0}'.format(
                            getattr(term, 'atom{0}'.format(n + 1)).index)
                    entry += '\n'
                    self._fh.write(entry)
                self._fh.write("</{0}>\n".format(directive[:-1]))

        # TODO: optional things
        self._fh.write("</configuration>\n")
        self._fh.write("</hoomd_xml>\n")

    def close(self):
        """Close the HOOMD xml file. """
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
    ethane.save("ethane.hoomdxml", forcefield='OPLS-AA')

