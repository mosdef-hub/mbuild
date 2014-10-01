import subprocess
import os


# TODO: could add a formatregister here for forcefield parameterizers
def topolbuild(traj, forcefield, charge=False, mol2_name='traj.mol2', top_name='traj.top'):
    """Use topolbuild to apply the OPLS-aa forcefield to a .pdb file.

    Args:
        traj (Trajectory): The trajectory to parameterize
    """
    traj.save(filename=mol2_name)
    mol2_name_base, _ = os.path.splitext(mol2_name)

    executable = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'topolbuild/src/topolbuild')
    ff_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'topolbuild/dat')

    options = [executable]
    # TODO: paths for other forcefields
    if forcefield == 'oplsaa':
        options.append('-dir {0}'.format(os.path.join(ff_directory, 'gromacs')))
        options.append('-ff oplsaa')

    if charge:
        assert forcefield == 'oplsaa'
        options.append('-charge')

    options.append('-n {0}'.format(mol2_name_base))

    # TODO: support for other options

    #pdb.set_trace()
    subprocess.call(options)

    #os.rename('ffoplsaa_TopolGen_.top', top_name)



