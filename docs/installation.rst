============
Installation
============

Install with conda (preferred)
------------------------------
::

    $ conda install --channel imodels mbuild


Install with pip
----------------
::

    $ pip install mbuild


Install from source
-------------------
::

    $ git clone https://github.com/iModels/mbuild
    $ cd mbuild
    $ python setup.py install

Dependencies
------------
To use mbuild, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit OS X Yosemite and Windows 7 machines.
        TravisCI is currently only set up to perform testing on Ubuntu 12.04
        LTS Server Edition 64 bit

    `Python <http://python.org>`_ = 2.7 or 3.3+
        TravisCI currently tests on 2.7, 3.3 and 3.4.

    `MDTraj <http://mdtraj.org/>`_ >=1.0.0
        MDTraj is a Python library for reading, writing and analyizing
        molecular dynamics trajectories. mBuild uses MDTraj as an entry and
        exit point for several types of molecule data formats. See their
        installation instructions
        `here <http://mdtraj.org/latest/installation.html>`_.

To make your life easier, we recommend that you use a pre-packaged Python
distribution like `Continuum's Anaconda <https://store.continuum.io/>`_
in order to get all of the dependencies.

Testing your installation
-------------------------

mBuild uses ``py.test`` for unit testing. To run them simply type run the
following while in the base directory::

    $ pip install pytest
    $ py.test

