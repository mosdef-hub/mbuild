#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]];   then MINICONDA=Miniconda3-latest-MacOSX-x86_64.sh; fi
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then MINICONDA=Miniconda3-latest-Linux-x86_64.sh;  fi

MINICONDA_MD5=$(curl -s https://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget https://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b
rm -f $MINICONDA

export PATH=$HOME/miniconda3/bin:$PATH
conda config --add channels omnia
conda config --add channels janschulz

conda update -yq conda

conda create -y -n myenv python=$PYTHON_VERSION
source activate myenv

conda install -y numpy
conda install -y scipy
conda install -y packmol 1.0.0 4
conda install -y nglview
conda install -y oset
conda install -y parmed
conda install -y mdtraj

conda install -y pytest
conda install -y jupyter
conda install -y nbformat
conda install -y ipykernel
conda install -y ipyext
