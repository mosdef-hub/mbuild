sudo apt-get update

MINICONDA=Miniconda-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b


export PATH=$HOME/miniconda/bin:$PAT
conda config --set always_yes yes --set changeps1 noH
conda install conda-build jinja2 binstar pip
conda config --add channels http://conda.binstar.org/omnia

conda install future coverage pytest jinja2
pip install coveralls
 pip install -e git://github.com/shirtsgroup/InterMol.git@forcetype#egg=intermol
python setup.py install