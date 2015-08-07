# Create the docs and push them to github pages
# ---------------------------------------------
conda install --yes sphinx numpydoc

python setup.py develop

cd docs/tutorials
ipython nbconvert --to html *.ipynb
cd ..
make html

source update_gh_pages.sh
