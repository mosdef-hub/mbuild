
def save_gromacs(traj, step=-1, top_file='mbuild.top', gro_file='mbuild.gro'):
    """Output a Trajectory as a GROMACS .gro and .top file.

    Args:
        traj:
        step:
        top_file:
        grofile:
    """

    with open(top_file, 'w') as f:


    with open(gro_file, 'w') as f:
        f.write('mbuild\n')
        f.write('{}\n'.format(traj.n_atoms))
        for n, data in enumerate(zip(traj.top.atoms, traj.xyz[step])):
            atom, xyz = data
            f.write('{:5d}{:<4s}{:6s}{:5d}'
                    '{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}\n'.format(
                atom.residue.index, atom.residue.name, atom.name, n + 1,
                xyz[0], xyz[1], xyz[2], 0.0, 0.0, 0.0))
        box = traj.unitcell_vectors[step]
        f.write('{:10.5f}{:10.5f}{:10.5f}'
                '{:10.5f}{:10.5f}{:10.5f}'
                '{:10.5f}{:10.5f}{:10.5f}'.format(
            box[0, 0], box[1, 1], box[2, 2],
            box[1, 0], box[2, 0], box[0, 1],
            box[2, 1], box[0, 2], box[1, 2]))
        f.write('\n')

if __name__ == "__main__":
    from mbuild.examples.ethane.ethane import Ethane
    ethane = Ethane().to_trajectory()
    save_gromacs(ethane)

