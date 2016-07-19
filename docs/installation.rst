============
Installation
============

Install with `conda <http://continuum.io/downloads>`_
-----------------------------------------------------
::

    $ conda install -c omnia -c janschulz -c imodels mbuild

.. note::
    The ``omnia`` channel provides `MDTraj <http://mdtraj.org/>`_,
    `ParmEd <http://parmed.github.io/ParmEd/>`_ and
    `packmol <http://www.ime.unicamp.br/~martinez/packmol/home.shtml/>`_,
    ``janschulz`` provides `ipyext <https://pypi.python.org/pypi/ipyext/0.1.0/>`_,
    and ``imodels`` is our channel where mBuild is hosted.

Alternatively you can add all the required channels to your ``.condarc``
after which you can simply ``conda install mbuild``::

    $ conda config --add channels omnia
    $ conda config --add channels janschulz
    $ conda config --add channels imodels
    $ conda install mbuild

.. note::
    The `MDTraj website <http://mdtraj.org/latest/new_to_python.html>`_ makes a
    nice case for using Python and in particular the
    `Anaconda scientific python distribution <http://continuum.io/downloads>`_
    to manage your numerical and scientific Python packages.

Install from source
-------------------
::

    $ git clone https://github.com/imodels/mbuild
    $ cd mbuild
    $ pip install -e .

To make your life easier, we recommend that you use a pre-packaged Python
distribution like `Continuum's Anaconda <https://store.continuum.io/>`_
in order to get all of the dependencies.

Testing your installation
-------------------------

mBuild uses ``py.test`` for unit testing. To run them simply type run the
following while in the base directory::

    $ conda install pytest
    $ py.test -v

