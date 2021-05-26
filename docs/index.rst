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
-------

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


Example system
---------------

Components in dashed boxes are drawn by hand using, e.g.,
`Avogadro <http://avogadro.cc>`_ or generated elsewhere. Each
component is wrapped as a simple python class with user defined attachment
sites, or ports. That's the hard part! Now mBuild can do the rest. Each component
further down the hierarchy is, again, a simple python class that describes
which piece should connect to which piece.

Ultimately, complex systems structures can be created with just a line or two
of code. Additionally, this approach seamlessly exposes tunable parameters within
the hierarchy so you can actually create whole families of structures simply
by adjusting a variable::

    pattern = Random2DPattern(20)  # A random arrangement of 20 pieces on a 2D surface.
    brush_layer = BrushLayer(chain_lenth=20, pattern=pattern, tile_x=3, tile_y=2)

.. image:: images/pmpc.png
    :align: center
    :alt: Zwitterionic brushes on beta-cristobalite substrate



Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.


.. toctree::
    	installation
    	docker
    	quick_start
    	quick_start/quick_start
    	structure_building_options
    	tutorials/tutorials
    	data_structures
	load_data
    	units
    	coordinate_transforms
    	recipes
	sim_engine_writers/sim_engine_writers
    	citing_mbuild
    	older_documentation
