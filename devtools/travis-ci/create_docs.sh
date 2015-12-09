# Create the docs and push them to github pages
# ---------------------------------------------
conda install --yes sphinx numpydoc mdtraj ipython-notebook oset networkx packmol
conda install --yes imolecule

python setup.py develop

cd docs/tutorials
ipython nbconvert --to html *.ipynb
cd ..
make html

source update_gh_pages.sh
