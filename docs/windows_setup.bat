# create a conda environment called mbuild
conda create -n iModels python=3 anaconda

# activate the environmnet
# windows
activate iModels
# UNIX
# . activate iModels

# add conda channels
conda config --add channels omnia
conda config --add channels patrickfuller

# install dependencies
conda install packmol nglview oset mdtraj pytest jupyter nbformat ipykernel imolecule

# clone mbuild
echo "****** git clone from source ******"
git clone https://github.com/iModels/mbuild.git
git clone https://github.com/iModels/foyer.git
git clone https://github.com/iModels/metamds.git
git clone https://github.com/ParmEd/ParmEd.git
git clone https://github.com/arose/nglview.git

echo "****** install from github source ******"
# cd to the mbuild directory
cd mbuild
# pip install it with the -e flag
pip install -e .

# install foyer
cd ..
cd foyer
pip install -e .

# install metamds
cd ..
cd metamds
pip install -e .

# install ParmEd
cd ..
cd ParmEd
pip install -e .

# install nglview
cd ..
cd nglview
pip install -e .
