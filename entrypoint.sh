#!/bin/sh

. /opt/conda/etc/profile.d/conda.sh
conda activate base
source activate mbuild-docker

if [ "$@" == "jupyter" ]; then
	jupyter notebook --no-browser --notebook-dir /home/anaconda/mbuild-notebooks --ip="0.0.0.0"
else
	$@
fi
