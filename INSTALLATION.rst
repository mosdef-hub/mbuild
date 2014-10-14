============
Installation
============

Install with pip
----------------

 $ pip install mbuild

Dependencies
------------
To use mbuild, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit Mac and Windows machines. TravisCI is
        currently only set up to perform testing on Debian.

    `Python <http://python.org>`_ == 2.7
         Once our unit tests flesh out a bit more, we intend to add support
         for >=2.6.

    `NumPy <http://numpy.scipy.org/>`_ >= 1.6.0
        Numpy is the base package for numerical computing in python.

    `MDTraj <http://mdtraj.org/>`_ >=1.0.0
        MDTraj is a Python library for reading, writing and analyizing
        molecular dynamics trajectories. mBuild uses MDTraj as an entry and
        exit point for all data.

We recommend that you use a pre-packaged Python distribution like
`Enthought's Canopy <https://www.enthought.com/products/canopy/>`_ or
`Continuum's Anaconda <https://store.continuum.io/>`_ in order to get all
of the dependencies.

=========================
Testing your installation
=========================