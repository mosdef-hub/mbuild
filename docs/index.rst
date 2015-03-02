.. mbuild documentation master file, created by
   sphinx-quickstart on Wed Oct  1 09:05:14 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

mBuild
======
*A hierarchical, component based molecule builder*

With just a few lines of mBuild code, you can assemble reusable components into
complex molecular systems for molecular dynamics simulations.


* mBuild is designed to minimize or even eliminate the need to explicitly translate and
  orient components when building systems: you simply tell it to connect two
  pieces!
* mBuild keeps track of the system's topology so you don't have to
  worry about keeping track of bonds when constructing chemically bonded
  structures from smaller components.

Example system
--------------

Components in dashed boxes are drawn by hand using, e.g.,
`Avogadro <http://avogadro.cc/wiki/Main_Page>`_ or generated elsewhere. Each
component is wrapped as a simple python class with user defined attachment
sites, or ports. That's the hard part! Now mBuild can do the rest. Each component
further down the hierarchy is, again, a simple python class that describes
which piece should connect to which piece.

Ultimately, complex systems structures can be created with just a line or two
of code. Additionally, this approach easily exposes tunable parameters within
the hierarchy so you can actually create whole families of structures simply
by adjusting a variable::

    mask = random_mask_2d(20)  # A random arrangement of 20 pieces on a 2D surface.
    brush_layer = BrushLayer(chain_lenth=20, mask=mask, tile_x=3, tile_y=2)

.. image:: images/pmpc.png
    :align: center
    :scale: 50%
    :alt: Zwitterionic brushes on beta-cristobalite substrate


.. image:: http://img.shields.io/badge/license-LGPLv2.1-blue.svg
    :target: https://www.gnu.org/licenses/lgpl.html

GNU LGPL version 2.1, or at your option a later version of the license.
Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.

Installation
------------
.. toctree::
   installation

Tutorials
---------
.. toctree::

   tutorials/tutorials

API Reference
-------------
.. toctree::
   :maxdepth: 1

   apidoc/mbuild
   apidoc/mbuild.tools
   apidoc/mbuild.formats

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
