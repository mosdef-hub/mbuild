sudo apt-get update
sudo apt-get install curl git libgmp-dev

MINICONDA=Miniconda-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b

export PATH=$HOME/miniconda/bin:$PATH
conda install --yes conda-build jinja2 binstar pip
conda config --add channels http://conda.anaconda.org/omnia
conda install --yes mdtraj

# rst -> ipynb -> html conversion
# pandoc (rst -> md)
conda install --yes --channel https://conda.anaconda.org/jsw-fnal pandoc
# notedown (md -> ipynb)
conda install --yes --channel https://conda.anaconda.org/auto notedown
# jsonschema: nbconvert and ipython requirements
conda install --yes jsonschema mistune

# libgmp hack
if [ ! -f /usr/lib/libgmp.so.3 ]; then
	sudo ln -s $(dpkg -L libgmp-dev | grep libgmp.so) /usr/lib/libgmp.so.3
fi