-----------------------
Quick Start
-----------------------

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

The `MoSDeF <https://mosdef.org>`_ software is comprised the following packages:
    * `mBuild <https://mbuild.mosdef.org/en/stable/>`_ -- A hierarchical, component based molecule builder
    * `foyer <https://foyer.mosdef.org/en/stable/>`_ -- A package for atom-typing as well as applying and disseminating forcefields
    * `GMSO <https://gmso.mosdef.org/en/stable/>`_ -- Flexible storage of chemical topology for molecular simulation

.. note::
    **foyer** and **GMSO** are used together with **mBuild** to create all the required files to conduct the simulations. Run time parameters for a simulation engine need to be created by the user.

In the following examples, different types of simulation boxes are constructed using the **MoSDeF** software.


Molecular simulations are usually comprised of many molecules contained in a
box (NPT and NVT ensembles), or boxes (GEMC and GCMC ensembles).
The **mBuild** library allows for easy generation of the simulation
box or boxes utilizing only a few lines of python code.


The following tutorials are available either as html or interactive `jupyter <https://jupyter.org/>`_ notebooks.


.. toctree::

   	load_files
   	Box_example
	fill_box_example
   	polymer_example
