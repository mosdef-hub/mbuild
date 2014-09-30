
def save_gromacs(traj, step=-1, basename='mbuild', forcefield='opls-aa'):
    """Output a Trajectory as a GROMACS .gro and .top file.

    Args:
        traj:
        step:
        top_file:
        grofile:
    """

    with open(basename + '.top', 'w') as f:
        if forcefield == 'opls-aa':
            f.write('#include "oplsaa.ff"\n\n')

        # TODO: iteration over chains
        f.write('\n[ moleculetype ]\n')
        f.write('; name nrexcl\n')
        f.write('mbuild 3\n')

        f.write('\n[ atoms ]\n')
        for atom in traj.topology.atoms:
            f.write('{:d} {:s} {:d} {:s} {:s} {:d} {:8.4f} {:8.4f}\n'.format(
                atom.index, atom.name, atom.residue.index, atom.residue.name,
                atom.name, 1, 0.0, atom.element.mass))

        f.write('\n[ bonds ]\n')
        for n, bond in enumerate(traj.topology.ff_bonds):
            f.write('{:d} {:d} {:d}\n'.format(
                n + 1, bond.atom1.index, bond.atom2.index))

        f.write('\n[ angles ]\n')
        for n, angle in enumerate(traj.topology.ff_angles):
            f.write('{:d} {:d} {:d} {:d}\n'.format(
                n + 1, angle.atom1.index, angle.atom2.index, angle.atom3.index))

        f.write('\n[ dihedrals ]\n')
        for n, dihedral in enumerate(traj.topology.ff_dihedrals):
            f.write('{:d} {:d} {:d} {:d}\n'.format(
                n + 1, dihedral.atom1.index, dihedral.atom2.index,
                dihedral.atom3.index, dihedral.atom4.index))

        f.write('\n[ system ]\n')
        f.write('{}\n'.format(basename))

        f.write('\n[ molecules ]\n')
        f.write('mbuild 1\n')

    with open(basename + '.gro', 'w') as f:
        f.write('{}\n'.format(basename))
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
    ethane.topology.find_forcefield_terms()
    save_gromacs(ethane, basename='ethane')

