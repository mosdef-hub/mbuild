============
Installation
============

Install with `conda <http://continuum.io/downloads>`_
-----------------------------------------------------
::

    $ conda install -c omnia -c mosdef -c conda-forge mbuild

Alternatively you can add all the required channels to your ``.condarc``
after which you can simply install without specifying the channels::

    $ conda config --add channels omnia
    $ conda config --add channels mosdef
    $ conda config --add channels conda-forge
    $ conda install mbuild

.. note::
    The `MDTraj website <http://mdtraj.org/latest/new_to_python.html>`_ makes a
    nice case for using Python and in particular the
    `Anaconda scientific python distribution <http://continuum.io/downloads>`_
    to manage your numerical and scientific Python packages.

Install with `pip <https://pypi.org/project/pip/>`_
---------------------------------------------------
::

    $ pip install mbuild

.. note::
    `PACKMOL <http://m3g.iqm.unicamp.br/packmol/>`_ is not available on pip
    but can be installed from source or via conda.

Install an editable version from source
---------------------------------------
::

    $ git clone https://github.com/mosdef-hub/mbuild
    $ cd mbuild
    $ pip install -e .

To make your life easier, we recommend that you use a pre-packaged Python
distribution like `Continuum's Anaconda <https://store.continuum.io/>`_
in order to get all of the dependencies.

Supported Python Versions
-------------------------

Python 3.6 and 3.7 are officially supported, including testing during
development and packaging. Support for Python 2.7 has been dropped as of 
August 6, 2019. Other Python versions, such as 3.8 and 3.5 and older, may
successfully build and function but no guarantee is made.

Testing your installation
-------------------------

mBuild uses ``py.test`` for unit testing. To run them simply type run the
following while in the base directory::

    $ conda install pytest
    $ py.test -v

