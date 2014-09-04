import pdb
from warnings import warn

from mdtraj.formats.registry import _FormatRegistry
from mbuild.compound import Compound
from mbuild.coordinate_transform import rotate_around_x

__all__ = ['load_hoomxml', 'save_hoomdxml']

@_FormatRegistry.register_loader('.hoomdxml')
def load_hoomdxml(filename, optional_nodes=False):
    """Load a HOOMD-blue XML file form disk.

    Args:
        filename (str): Path to xml file.
        optional_nodes(bool, optional):

    Returns:
        traj (md.Trajectory):
    """

    import xml.etree.cElementTree as etree

    import numpy as np
    import pandas as pd

    #from mbuild.trajectory import Trajectory
    from mdtraj.core.trajectory import Trajectory
    #from mbuild.topology import Topology
    from mdtraj.core.topology import Topology

    tree = etree.parse(filename)

    config = tree.getroot().find('configuration')
    # Required nodes for valid HOOMD simulation: box, position and type.
    box = config.find('box')
    lx = float(box.attrib['lx'])
    ly = float(box.attrib['ly'])
    lz = float(box.attrib['lz'])
    xy = float(box.attrib['xy'])
    xz = float(box.attrib['xz'])
    yz = float(box.attrib['yz'])
    unitcell_vectors = np.array([[[lx,  xy*ly, xz*lz],
                                  [0.0, ly,    xy*lz],
                                  [0.0, 0.0,   lz   ]]])

    xyz = list()
    for n, pos in enumerate(config.find('position').text.splitlines()[1:]):
        xyz.append([float(x) for x in pos.split()])
    n_atoms = n + 1

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
                    parsed_line = [int(x) if x.isdigit() else float(x) for x in raw_line.split()]
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

    # found = [node for node in optional_data.keys() if node != 'per_particle']
    # found.extend([node for node in per_particle_df.columns])
    # print "Parsed the following optional nodes from '{0}':".format(filename)
    # print found

    atoms_df = pd.DataFrame(atom_types, columns=['name'])
    atoms_df['element'] = ['' for n in range(n_atoms)]
    atoms_df['resSeq'] = np.ones(n_atoms, dtype='int')
    atoms_df['resName'] = ['RES' for n in range(n_atoms)]
    atoms_df['chainID'] = np.ones(n_atoms, dtype='int')

    # TODO: Infer chains by finding isolated bonded structures.
    if 'bond' in optional_data:
        bonds = optional_data['bond'][['id0', 'id1']].values
    else:
        bonds = np.empty(shape=(0,2), dtype="int")
    top = Topology.from_dataframe(atoms_df, bonds=bonds)

    traj = Trajectory(xyz=np.array(xyz, dtype=np.float64), topology=top)
    traj.unitcell_vectors = unitcell_vectors
    if optional_nodes:
        return traj, optional_data
    return traj

# TODO: decide if we want this to be a class
def save_hoomdxml(traj, step=-1, optional_nodes=None, filename='mbuild.xml'):
    """Output a Compound as a HOOMD XML file.

    Args:
        component (Compound): The Compound to be output.
        filename (str, optional): Path of the output file.

    """
    with open(filename, 'w') as xml_file:
        xml_file.write("""<?xml version="1.3" encoding="UTF-8"?>\n""")
        xml_file.write("""<hoomd_xml>\n""")
        xml_file.write("""<configuration time_step="0">\n""")

        lx, ly, lz = traj.unitcell_lengths[0]
        xy = traj.unitcell_vectors[0,1,0] / ly
        xz = traj.unitcell_vectors[0,2,0] / lz
        yz = traj.unitcell_vectors[0,2,1] / lz

        xml_file.write("""<box lx="{0}" ly="{1}" lz="{2}" xy="{3}" xz="{4}" yz="{5}" />\n""".format(lx, ly, lz, xy, xz, yz))

        xml_file.write("""<position num="{0}">\n""".format(traj.n_atoms))
        for atom in traj.xyz[step]:
            x, y, z = atom
            xml_file.write("{0:8.5f} {1:8.5f} {2:8.5f}\n".format(x, y, z))
        xml_file.write("</position>\n")

        xml_file.write("""<type num="{0}">\n""".format(traj.n_atoms))
        for atom in traj.top.atoms:
            xml_file.write("{0}\n".format(atom.name))
        xml_file.write("</type>\n")

        # TODO: optional things
        xml_file.write("</configuration>\n")
        xml_file.write("</hoomd_xml>\n")

if __name__ == "__main__":
    # traj, optional_data = load_hoomdxml('init.xml', optional_nodes=True)
    # save_hoomdxml(traj, filename='init_out.xml')
    # pdb.set_trace()


    from examples.ethane.ethane import Ethane
    ethane = Ethane()
    import numpy as np
    rotate_around_x(ethane, np.pi)
    ethane.save("ethane.hoomdxml")

    ethane.update_from_file("ethane.hoomdxml")

    c = Compound()
    c.append_from_file("ethane.hoomdxml")

    for atom in c.atoms():
        print atom