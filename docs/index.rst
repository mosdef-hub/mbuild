mBuild
=======

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

*A hierarchical, component based molecule builder*

With just a few lines of mBuild code, you can assemble reusable components into
complex molecular systems for molecular simulations.


* mBuild is designed to minimize or even eliminate the need to explicitly translate and
  orient components when building systems: you simply tell it to connect two
  pieces!
* mBuild keeps track of the system's topology so you don't have to
  worry about manually defining bonds when constructing chemically bonded
  structures from smaller components.



mBuild is a part of the MoSDeF ecosystem
----------------------------------------

The **mBuild** software, in conjunction with the other `Molecular Simulation Design Framework (MoSDeF) <https://mosdef.org>`_ tools, supports a wide range of
simulation engines, including `Cassandra <https://cassandra.nd.edu>`_, `GPU Optimized Monte Carlo (GOMC) <http://gomc.eng.wayne.edu>`_, `GROMACS <https://www.gromacs.org>`_,
`HOOMD-blue <http://glotzerlab.engin.umich.edu/hoomd-blue/>`_, and
`Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) <https://lammps.sandia.gov>`_.
The **mBuild** and **MoSDeF** tools allow simulation reproducibility
across the various simulation engines, eliminating the need to be an expert user in all
the engines to replicate, continue, or advance the existing research. Additionally,
the software can auto-generate many different systems,
allowing large-scale screening of chemicals and materials using
`Signac <https://signac.io>`_ to manage the simulations and data.

The `MoSDeF <https://mosdef.org>`_ software is comprised the following packages:
    * `mBuild <https://mbuild.mosdef.org/en/stable/>`_ -- A hierarchical, component based molecule builder
    * `foyer <https://foyer.mosdef.org/en/stable/>`_ -- A package for atom-typing as well as applying and disseminating forcefields
    * `GMSO <https://gmso.mosdef.org/en/stable/>`_ -- Flexible storage of chemical topology for molecular simulation


.. toctree::
	:caption: Getting Started
    	:maxdepth: 2

	getting_started/example_system
	getting_started/installation/installation_toc
    	getting_started/quick_start/quick_start
	getting_started/writers/writers
	getting_started/tutorials/tutorials

.. toctree::
	:caption: Topic Guides
    	:maxdepth: 2

	topic_guides/recipe_development
    	topic_guides/data_structures
	topic_guides/load_data
    	topic_guides/coordinate_transforms
    	topic_guides/recipes

.. toctree::
    	:caption: Reference
    	:maxdepth: 2

	reference/units
    	reference/citing_mbuild
    	reference/older_documentation

.. toctree::
    :caption: API Documentation

    autoapi/mbuild/index
