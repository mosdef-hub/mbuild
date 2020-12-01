## mBuild: a hierarchical, component based molecule builder

[![Gitter chat](https://badges.gitter.im/mosdef-hub/gitter.svg)](https://gitter.im/mosdef-hub/Lobby)
[![AZP Build Status](https://dev.azure.com/mosdef/mbuild/_apis/build/status/mosdef-hub.mbuild?branchName=master)](https://dev.azure.com/mosdef/mbuild/_build/latest?definitionId=1&branchName=master)
[![PyPI Version](https://badge.fury.io/py/mbuild.svg)](https://pypi.python.org/pypi/mbuild)
[![Anaconda Badge](https://anaconda.org/conda-forge/mbuild/badges/version.svg)](https://anaconda.org/conda-forge/mbuild)
[![codecov](https://codecov.io/gh/mosdef-hub/mbuild/branch/master/graph/badge.svg)](https://codecov.io/gh/mosdef-hub/mbuild)

With just a few lines of mBuild code, you can assemble reusable components into
complex molecular systems for molecular dynamics simulations.

* mBuild is designed to minimize or even eliminate the need to explicitly translate and
  orient components when building systems: you simply tell it to connect two
  pieces!
* mBuild keeps track of the system's topology so you don't have to
  worry about manually defining bonds when constructing chemically bonded
  structures from smaller components.

To learn more, get started or contribute, check out our [website](http://mbuild.mosdef.org).

If you use this package, please cite [our paper](http://dx.doi.org/10.1007/978-981-10-1128-3_5
). The BibTeX reference is
```
@article{Klein2016mBuild,
      author = "Klein, Christoph and Sallai, János and Jones, Trevor J. and Iacovella, Christopher R. and McCabe, Clare and Cummings, Peter T.",
      title = "A Hierarchical, Component Based Approach to Screening Properties of Soft Matter",
      booktitle = "Foundations of Molecular Modeling and Simulation",
      series = "Molecular Modeling and Simulation: Applications and Perspectives",
      year = "2016",
      doi = "http://dx.doi.org/10.1007/978-981-10-1128-3_5" 
}
```

#### Quick Start with Docker
To use `mbuild` in a jupyter-notebook that runs from a docker container with all the dependencies installed use the following command:

```sh
$ docker pull mosdef/mbuild:latest
$ docker run -it --name mbuild -p 8888:8888 mosdef/mbuild:latest su anaconda -s\
      /bin/sh -l -c "jupyter-notebook --no-browser --ip="0.0.0.0" --notebook-dir\
      /home/anaconda/mbuild-notebooks
```

Alternatively, you can also start a Bourne shell directly:
```sh
$ docker run -it --name mbuild mosdef/mbuild:latest
```

To learn more about using `mBuild` with docker, please refer to the documentation [here](https://mbuild.mosdef.org/en/latest/docker.html).

#### Tutorials

*Interactive tutorials can be found here:* 

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/mosdef-hub/mbuild_tutorials/master)

Components in dashed boxes are drawn by hand using, e.g.,
[Avogadro](https://avogadro.cc/) or generated elsewhere. Each
component is wrapped as a simple python class with user defined attachment
sites, or ports. That's the "hard" part! Now mBuild can do the rest. Each component
further down the hierarchy is, again, a simple python class that describes
which piece should connect to which piece.

Ultimately, complex structures can be created with just a line or two
of code. Additionally, this approach seamlessly exposes tunable parameters within
the hierarchy so you can actually create whole families of structures simply
by adjusting a variable:

```python
import mbuild as mb
from mbuild.examples import PMPCLayer

pattern = mb.Random2DPattern(20)  # A random arrangement of 20 pieces on a 2D surface.
pmpc_layer = PMPCLayer(chain_length=20, pattern=pattern, tile_x=3, tile_y=2)
```

![Zwitterionic brushes on beta-cristobalite substrate](docs/images/pmpc.png)
#### Community Recipes
Use case-specific systems can be generated via mBuild recipes.
Some users have graciously contributed recipes for particular systems, including:

* [Graphene slit pores](https://github.com/rmatsum836/Pore-Builder)
* [Nanodroplets on graphene](https://github.com/ftiet/droplet-builder)
* [Coarse-grained DNA](https://github.com/zijiewu3/mbuild_ONA)
* [Lipid bilayers](https://github.com/uppittu11/mbuild_bilayer)

#### [![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.

This material is based upon work supported by the National Science Foundation under grants NSF CBET-1028374 and NSF ACI-1047828. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
