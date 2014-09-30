from mdtraj.utils import in_units_of


def save_xyz(traj, step=-1, filename='mbuild.xyz'):
    """ """
    with open(filename, 'w') as f:
        f.write(str(traj.n_atoms) + '\n\n')
        for atom, xyz in zip(traj.top.atoms, traj.xyz[step]):
            x, y, z = in_units_of(xyz, 'nanometers', 'angstroms')
            f.write("{0} {1} {2} {3}\n".format(atom.name, x, y, z))
