import numpy as np
import pandas as pd
from xml.etree import cElementTree
from mdtraj.formats.registry import _FormatRegistry

from mbuild.compound import Compound


__all__ = ['load_hoomdxml', 'save_hoomdxml']


@_FormatRegistry.register_loader('.hoomdxml')
def load_hoomdxml(filename, optional_nodes=True, lj_units=None):
    """Load a HOOMD-blue XML file form disk.

    Note: lj_units need to be normalized by nm, kJ/mol, and amu

    Required nodes for valid HOOMD simulation: box, position and type.

    Parameters
    ----------
    filename : str
        Path to xml file.
    optional_nodes : bool, optional
        Read all nodes in file other than 'box', 'position' and 'type'.

    Returns
    -------
    compound : mb.Compound

    """
    # Get fundamental LJ units.
    if lj_units is None:
        lj_units = {'distance': 1.0,
                    'energy': 1.0,
                    'mass': 1.0}
    else:
        assert isinstance(lj_units, dict)
        assert 'distance' in lj_units
        assert 'energy' in lj_units
        assert 'mass' in lj_units

    # Other derived LJ units.
    lj_units['time'] = (np.sqrt(lj_units['mass'] * lj_units['distance']**2.0
                        / lj_units['energy']))
    lj_units['velocity'] = lj_units['distance'] / lj_units['time']
    lj_units['acceleration'] = lj_units['distance'] / lj_units['time']**2.0
    lj_units['diameter'] = lj_units['distance']
    lj_units['charge'] = 1.0
    # TODO: figure out charge
    lj_units['moment_inertia'] = lj_units['mass'] * lj_units['distance']**2.0
    lj_units['image'] = 1.0
    lj_units['body'] = 1.0
    lj_units['orientation'] = 1.0

    compound = Compound()
    tree = cElementTree.parse(filename)
    config = tree.getroot().find('configuration')

    # Unitcell info.
    box = config.find('box')
    unitcell_vectors = _unitcell_vectors(box) * lj_units['distance']

    # Coordinates.
    xyz = list()
    for pos in config.find('position').text.splitlines()[1:]:
        xyz.append([float(x) * lj_units['distance'] for x in pos.split()])
    n_atoms = len(xyz)

    # TODO: Create custom elements based on type names. Probably want a prefix
    #       or suffix to avoid overwriting atomic elements.
    atom_types = list()
    for atom_type in config.find('type').text.splitlines()[1:]:
        atom_types.append(atom_type)

    optional_data = dict()
    if optional_nodes:
        # Create a dataframe with all available per-particle information.
        per_particle_df = pd.DataFrame()
        per_particle_nodes = ['image', 'velocity', 'acceleration', 'mass',
                'diameter', 'charge', 'body', 'orientation', 'moment_inertia']
        for node in per_particle_nodes:
            parsed_node_text = list()
            try:
                node_text = config.find(node).text.splitlines()[1:]
                for raw_line in node_text:
                    # TODO: not robust when e.g. charges are provided as ints
                    parsed_line = [int(x) if x.isdigit() else float(x) * lj_units[node] 
                            for x in raw_line.split()]
                    if len(parsed_line) == 1:
                        parsed_line = parsed_line[0]
                    parsed_node_text.append(parsed_line)
                per_particle_df[node] = parsed_node_text
            except AttributeError as err:
                pass
        if not per_particle_df.empty:
            optional_data['per_particle'] = per_particle_df

        # Add a dataframe for each available multi-particle node.
        multi_particle_nodes = [('bond', 2), ('angle', 3), ('dihedral', 4),
                                ('improper', 4)]
        for node, n_indices in multi_particle_nodes:
            parsed_node_text = list()
            try:
                node_text = config.find(node).text.splitlines()[1:]
                for raw_line in node_text:
                    parsed_line = [int(x) if x.isdigit() else x for x in raw_line.split()]
                    parsed_node_text.append(parsed_line)
                columns = ['id'+str(n) for n in range(n_indices)]
                columns.insert(0, node+'type')
                multi_particle_df = pd.DataFrame(parsed_node_text, columns=columns)
                optional_data[node] = multi_particle_df
            except (AttributeError, ValueError) as err:
                pass

        # TODO: read wall
        # <wall> has its information stored as attributes

    atoms_df = pd.DataFrame(atom_types, columns=['name'])
    atoms_df['serial'] = range(n_atoms)
    atoms_df['element'] = ['' for _ in range(n_atoms)]
    atoms_df['resSeq'] = np.ones(n_atoms, dtype='int')
    atoms_df['resName'] = ['RES' for _ in range(n_atoms)]
    atoms_df['chainID'] = np.ones(n_atoms, dtype='int')
    # TODO: Infer chains by finding isolated bonded structures.
    if 'bond' in optional_data:
        bonds = optional_data['bond'][['id0', 'id1']].values
    else:
        bonds = np.empty(shape=(0, 2), dtype="int")

        
    return compound


def _unitcell_vectors(box):
    """Parse unitcell vectors from box node. """

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

    return np.array([[[lx,   xy*ly, xz*lz],
                       [0.0,    ly, yz*lz],
                       [0.0,   0.0,    lz]]])


def save_hoomdxml(traj, step=-1, optional_nodes=None, filename='mbuild.xml'):
    """Output a Trajectory as a HOOMD XML file.

    Args:
        traj (md.Trajectory): The Trajectory to be output.
        filename (str, optional): Path of the output file.

    """
    with open(filename, 'w') as xml_file:
        xml_file.write("""<?xml version="1.3" encoding="UTF-8"?>\n""")
        xml_file.write("""<hoomd_xml>\n""")
        xml_file.write("""<configuration time_step="0">\n""")

        lx, ly, lz = traj.unitcell_lengths[step]
        xy = traj.unitcell_vectors[0, 1, 0] / ly
        xz = traj.unitcell_vectors[0, 2, 0] / lz
        yz = traj.unitcell_vectors[0, 2, 1] / lz

        xml_file.write("""<box lx="{0}" ly="{1}" lz="{2}" xy="{3}" xz="{4}" yz="{5}" />\n""".format(lx, ly, lz, xy, xz, yz))

        xml_file.write("""<position num="{0}">\n""".format(traj.n_atoms))
        for atom in traj.xyz[step]:
            x, y, z = atom
            xml_file.write("{0:8.5f} {1:8.5f} {2:8.5f}\n".format(float(x), float(y), float(z)))
        xml_file.write("</position>\n")

        xml_file.write("""<type num="{0}">\n""".format(traj.n_atoms))
        for atom in traj.top.atoms:
            xml_file.write("{0}\n".format(atom.atomtype))
        xml_file.write("</type>\n")

        optional_directives = [('bonds', 2), ('angles', 3), ('dihedrals', 4),
                               ('impropers', 4)]
        for directive, n_terms in optional_directives:
            if getattr(traj.top, '_ff_{0}'.format(directive)):
                xml_file.write("""<{0}>\n""".format(directive[:-1]))
                for term in getattr(traj.top, 'ff_{0}'.format(directive)):
                    entry = '{0} '.format(term.kind)
                    for n in range(n_terms):
                        entry += '{0} '.format(
                            getattr(term, 'atom{0}'.format(n + 1)).index)
                    entry += '\n'
                    xml_file.write(entry)
                xml_file.write("</{0}>\n".format(directive[:-1]))

        # TODO: optional things
        xml_file.write("</configuration>\n")
        xml_file.write("</hoomd_xml>\n")

if __name__ == "__main__":
    from mbuild.testing.tools import get_fn

    import numpy as np

    from mbuild.examples.ethane.ethane import Ethane
    from mbuild.coordinate_transform import rotate_around_x

    ethane = Ethane()
    rotate_around_x(ethane, np.pi)
    ethane = ethane.to_trajectory()
    ethane.top.find_forcefield_terms()
    ethane.save("ethane.hoomdxml")

    # ethane.update_from_file("ethane.hoomdxml")
    #
    # from mbuild.compound import Compound
    # c = Compound()
    # c.append_from_file("ethane.hoomdxml")
    #
    # for atom in c.atoms():
    #     print atom
