# create a conda environment called mbuild
conda create -n mbuild python=3 anaconda

# activate the environmnet
# windows
activate mbuild
# UNIX
#. activate mbuild

# add conda channels
conda config --add channels omnia

# install dependencies
conda install packmol nglview oset parmed mdtraj pytest jupyter nbformat ipykernel

# clone mbuild
git clone git@github.com:iModels/mbuild.git

# cd to the mbuild directory
cd mbuild
# pip install it with the -e flag
pip install -e .
