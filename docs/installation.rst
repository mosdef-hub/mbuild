============
Installation
============

Install with `conda <https://repo.anaconda.com/miniconda/>`_
-----------------------------------------------------
::

    $ conda install -c conda-forge -c omnia mbuild

Alternatively you can add all the required channels to your ``.condarc``
after which you can simply install without specifying the channels::

    $ conda config --add channels omnia
    $ conda config --add channels conda-forge
    $ conda install mbuild

.. note::
    The order in which channels are added matters: ``conda-forge`` should be the highest priority as a result of being added last. In your ``.condarc`` file, it should be listed first.

.. note::
    The `MDTraj website <http://mdtraj.org/1.9.3/new_to_python.html>`_ makes a
    nice case for using Python and in particular the
    `Anaconda scientific python distribution <https://www.anaconda.com/products/individual>`_
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
distribution like `Continuum's Anaconda <https://www.anaconda.com/products/individual/>`_
in order to get all of the dependencies.

Supported Python Versions
-------------------------

Python 3.6, 3.7 and 3.8 are officially supported, including testing during
development and packaging. Support for Python 2.7 has been dropped as of
August 6, 2019. Other Python versions, such as 3.9 and 3.5 and older, may
successfully build and function but no guarantee is made.

Testing your installation
-------------------------

mBuild uses ``py.test`` for unit testing. To run them simply type run the
following while in the base directory::

    $ conda install pytest
    $ py.test -v

