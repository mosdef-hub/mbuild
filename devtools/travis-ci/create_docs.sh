# Create the docs and push them to github pages
# ---------------------------------------------
conda install --yes sphinx numpydoc

python setup.py develop

cd docs
make html

source update_gh_pages.sh
