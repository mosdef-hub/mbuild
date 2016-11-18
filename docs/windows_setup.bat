@echo off

rem create a conda environment called mbuild
rem conda create -n iModels python=3 anaconda

rem activate the environmnet
rem windows
rem activate iModels
rem UNIX
rem . activate iModels

rem add conda channels
conda config --add channels omnia
conda config --add channels patrickfuller

rem install dependencies
conda install packmol nglview oset mdtraj pytest jupyter nbformat ipykernel imolecule

rem clone mbuild

text
****** git clone from source ******
endtext

git clone https://github.com/iModels/mbuild.git
git clone https://github.com/iModels/foyer.git
git clone https://github.com/iModels/metamds.git
git clone https://github.com/ParmEd/ParmEd.git
git clone https://github.com/arose/nglview.git

text
****** install from github source ******
endtext

rem cd to the mbuild directory
cd mbuild
rem pip install it with the -e flag
pip install -e .

rem install foyer
cd ..
cd foyer
pip install -e .

rem install metamds
cd ..
cd metamds
pip install -e .

rem install ParmEd
cd ..
cd ParmEd
pip install -e .

rem install nglview
cd ..
cd nglview
pip install -e .
