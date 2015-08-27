echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


if [[ "2.7 3.4" =~ "$python" ]]; then
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package mbuild-dev $HOME/miniconda/conda-bld/linux-64/mbuild-*
    conda convert $HOME/miniconda/conda-bld/linux-64/mbuild-* -p all
    ls
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package mbuild-dev linux-32/mbuild-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package mbuild-dev win-32/mbuild-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package mbuild-dev win-64/mbuild-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package mbuild-dev osx-64/mbuild-*
fi

if [[ "$python" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${python}"; exit 0
fi
