#! usr/bin/env bash

# This script is inspired from frued's pypi publish script at
# https://github.com/glotzerlab/freud/blob/v2.3.0/.circleci/deploy-linux.sh

cat << EOF > ~/.pypirc
[distutils]
index-servers=
    pypi
    testpypi
[pypi]
username: ${PYPI_USERNAME}
password: ${PYPI_PASSWORD}
[testpypi]
repository: https://test.pypi.org/legacy/
username: ${PYPI_USERNAME}
password: ${PYPI_PASSWORD}
EOF

source activate test-environment
pip install setuptools twine wheel
python setup.py sdist bdist_wheel


if [ -z $1 ]; then
    echo "A repository (\"pypi\" or \"testpypi\") must be provided as the first argument."
    exit 1
fi

twine upload --skip-existing -r $1 dist/*
