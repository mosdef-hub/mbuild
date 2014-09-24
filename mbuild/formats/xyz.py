__author__ = 'CTK'


def save_xyz(traj, step=-1, filename='mbuild.xyz'):
    """ """
    with open(filename, 'w') as f:
        f.write(str(traj.n_atoms) + '\n\n')
        for atom, xyz in zip(traj.top.atoms, traj.xyz[step]):
            # TODO: proper unit conversion
            x, y, z = xyz * 10.0
            f.write("{0} {1} {2} {3}\n".format(atom.name, x, y, z))
