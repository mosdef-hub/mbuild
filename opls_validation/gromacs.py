from mbuild.trajectory import Trajectory


def load_top(filename):
    pdb_file = "{}-gas.pdb".format(filename[:-4])
    traj = Trajectory.load(pdb_file)

    opls_types = list()

    with open(filename, 'r') as top:
        current = None
        for line in top:
            if ';' in line:
                line = line[:line.index(';')]
            stripped = line.strip()
            if stripped.startswith('*') or len(stripped) == 0:
                continue
            elif stripped.startswith('['):
                current = stripped[1:-1].strip()
            elif current == 'atoms':
                opls_types.append(stripped.split()[1])
            elif current == 'bonds':
                atom1 = traj.topology.atom(int(stripped.split()[0]) - 1)
                atom2 = traj.topology.atom(int(stripped.split()[1]) - 1)
                traj.topology.add_bond(atom1, atom2)
    compound = traj.to_compound()
    return compound, opls_types

