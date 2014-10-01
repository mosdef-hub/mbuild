import subprocess
import os

# TODO: could add a formatregister here for forcefield parameterizers
def topolgen(traj, pdb_name='traj.pdb', top_name='traj.top'):
    """Use topolgen to apply the OPLS-aa forcefield to a .pdb file.

    Args:
        traj (md.Trajectory): The trajectory to parameterize
        top_name (str, optional): Rename the output .top file to something
        other than the default 'ffoplsaa_TopolGen_.top'
    """
    traj.save(pdb_name)
    directory = os.path.dirname(os.path.realpath(__file__))
    perl_file = os.path.join(directory, 'topolgen/topolgen.pl')
    pipe = subprocess.Popen(["perl", perl_file, pdb_name])
    pipe.wait()

    os.rename('ffoplsaa_TopolGen_*.top', top_name)



