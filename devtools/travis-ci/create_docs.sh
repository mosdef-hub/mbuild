#!/bin/bash

# Create the docs
# get the directory in which this script is located
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"


# Print each line, exit on error
set -ev

# we assume that we're running in a virtualenv with all mbuild already installed and tested

# install additional packages required for doc generation
conda install --yes sphinx numpydoc sphinx sphinx_rtd_theme widgetsnbextension ipywidgets
pip install mock

pushd $DIR/../..

# this will generate mbuild/versions.py
python setup.py --name

cd docs

# clean leftovers from previous run
rm -f *.mol2
rm -f *.xyz
rm -f *.pdb

# FIXME: the docs folder is the working directory for the tutorial notebooks
# so we need to symlink the resources here from tutorials/ 

ln -s tutorials/ch3.pdb .

make html
