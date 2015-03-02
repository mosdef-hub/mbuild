============
Installation
============

Install with pip
----------------
::

    $ pip install mbuild

Dependencies
------------
To use mbuild, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit Mac and Windows machines. TravisCI is
        currently only set up to perform testing on Debian.

    `Python <http://python.org>`_ == 2.7
         Once our unit tests flesh out a bit more, we intend to add support
         for 3.4+

    `MDTraj <http://mdtraj.org/>`_ >=1.0.0
        MDTraj is a Python library for reading, writing and analyizing
        molecular dynamics trajectories. mBuild uses MDTraj as an entry and
        exit point for molecule data.

Optional packages:

    `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_
        VMD (Visual Molecular Dynamics) is a visualization program. Currently,
        we only use a (very) thin wrapper to call it from the command-line
        when you want to look at a Compound that you've built. A more robust
        integration with partial control from Python is planned in the near
        future.

To make your life easier, we recommend that you use a pre-packaged Python
distribution like `Enthought's Canopy <https://www.enthought.com/products/canopy/>`_
or `Continuum's Anaconda <https://store.continuum.io/>`_ in order to get all
of the dependencies.

Testing your installation
-------------------------

mBuild uses `py.test` for unit testing. To run them simply type run the
following while in the base directory::

    $ py.test

We need a LOT more tests so any help here is especially welcome!
