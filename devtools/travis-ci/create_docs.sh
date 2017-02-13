# Create the docs

# Print each line, exit on error
set -ev

# we assume that we're running in a virtualenv with all mbuild already installed and tested

# install additional packages required for doc generation
conda install --yes sphinx numpydoc sphinx sphinx_rtd_theme widgetsnbextension ipywidgets

pushd docs

# clean leftovers from previous run
rm -rf _build
rm -f *.mol2
rm -f *.xyz
rm -f *.pdb

# FIXME: the docs folder is the working directory for the tutorial notebooks
# so we need to symlink the resources here from tutorials/ 

ln -s tutorials/ch3.pdb .

make html

popd