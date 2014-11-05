import os

from pkg_resources import resource_filename

from mbuild.trajectory import Trajectory


def get_fn(name):
    """Get the full path to one of the reference files shipped for testing

    In the source distribution, these files are in ``mbuild/testing/reference``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Args:
        name (str): Name of the file to load (with respect to the reference/ folder).

    """
    fn = resource_filename('mbuild', os.path.join('testing', 'reference', name))

    if not os.path.exists(fn):
        raise ValueError('Sorry! %s does not exists. If you just '
                         'added it, you\'ll have to re install' % fn)

    return fn


def load_top_opls(toppath):
    """
    Args:
        toppath (str): Path to the top file.
    """
    split_path = os.path.split(toppath)
    filename = split_path[-1]
    pdb_file = "{}-gas.pdb".format(filename[:-4])
    pdb_path = os.path.join(split_path[0], pdb_file)
    traj = Trajectory.load(pdb_path)

    opls_types = list()
    with open(toppath, 'r') as top:
        current = None
        for line in top:
            if ';' in line:
                line = line[:line.index(';')]
            stripped = line.strip()
            if stripped.startswith('*') or len(stripped) == 0:
                continue
            elif stripped.startswith('['):
                current = stripped[1:-1].strip()
            elif current == 'moleculetype':
                mol_name = stripped.split()[0]
            elif current == 'atoms':
                opls_string = stripped.split()[1]
                if 'opls' not in opls_string:
                    print "Found non 'opls' type in {} ({}): {}.".format(
                        mol_name, filename, opls_string)
                    print "Ignoring file. Need to come up with proper support.\n"
                    return
                opls_types.append(opls_string.split('_')[1])
            elif current == 'bonds':
                atom1 = traj.topology.atom(int(stripped.split()[0]) - 1)
                atom2 = traj.topology.atom(int(stripped.split()[1]) - 1)
                traj.topology.add_bond(atom1, atom2)
    compound = traj.to_compound()
    return compound, opls_types, mol_name
