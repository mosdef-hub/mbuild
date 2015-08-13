============
Installation
============

Install with conda (preferred)
------------------------------
::

    $ conda install -c omnia -c patrickfuller -c imodels mbuild

or alternatively you can add all the required channels to your ``.condarc``
after which you can simply ``conda install mbuild``::

    $ conda config --add channels omnia
    $ conda config --add channels patrickfuller
    $ conda config --add channels imodels
    $ conda install mbuild

.. note::
    The ``omnia`` channel provides `MDTraj <http://mdtraj.org/>`_,
    ``patrickfuller`` provides `imolecule <http://patrick-fuller.com/imolecule/>`_,
    and ``imodels`` is our channel where mBuild is hosted.

Install from source
-------------------
::

    $ git clone https://github.com/imodels/mbuild
    $ cd mbuild
    $ python setup.py install

Dependencies
------------
To use mBuild, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit OS X and Windows machines.
        Continuous integration is run on TravisCI and Appveyor.

    `Python <http://python.org>`_ = 2.7 or 3.3+
        TravisCI and Appveyor currently test on 2.7 and 3.4.

    `MDTraj <http://mdtraj.org/>`_ >=1.0.0
        MDTraj is a Python library for reading, writing and analyizing
        molecular dynamics trajectories. mBuild uses MDTraj as an entry and
        exit point for several types of molecule data formats. See their
        installation instructions
        `here <http://mdtraj.org/latest/installation.html>`_.

    `imolecule <http://patrick-fuller.com/imolecule/>`_ >=0.1.9
        imolecule is an embeddable webGL molecule viewer and file format
        converter. mBuild uses imolecule to visualize components either inline
        if you're using an IPython notebook or in a browser window if not.


To make your life easier, we recommend that you use a pre-packaged Python
distribution like `Continuum's Anaconda <https://store.continuum.io/>`_
in order to get all of the dependencies.

Testing your installation
-------------------------

mBuild uses ``py.test`` for unit testing. To run them simply type run the
following while in the base directory::

    $ conda install pytest
    $ py.test -v

