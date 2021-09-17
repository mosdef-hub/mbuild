#!/bin/sh

. /opt/conda/etc/profile.d/conda.sh
conda activate base
conda activate mbuild-dev

if [ "$@" == "jupyter" ]; then
	jupyter notebook --no-browser --notebook-dir /home/anaconda/data --ip="0.0.0.0"
else
	$@
fi
