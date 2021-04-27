============
Installation
============

Install with `conda <https://repo.anaconda.com/miniconda/>`_
-----------------------------------------------------
::

    $ conda install -c conda-forge mbuild

Alternatively you can add all the required channels to your ``.condarc``
after which you can simply install without specifying the channels::

    $ conda config --add channels conda-forge
    $ conda install mbuild

.. note::
    The order in which channels are added matters: ``conda-forge`` should be the highest priority as a result of being added last. In your ``.condarc`` file, it should be listed first.

.. note::
    The `MDTraj website <http://mdtraj.org/1.9.3/new_to_python.html>`_ makes a
    nice case for using Python and in particular the
    `Anaconda scientific python distribution <https://www.anaconda.com/products/individual>`_
    to manage your numerical and scientific Python packages.

Install an editable version from source
---------------------------------------

To make your life easier, we recommend that you use a pre-packaged Python
distribution like `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
in order to get all of the dependencies::

    $ git clone https://github.com/mosdef-hub/mbuild
    $ cd mbuild
    $ conda env create -f environment-dev.yml
    $ conda activate mbuild-dev
    $ pip install -e .

.. note::
    The above installation is for OSX and Unix. If you are using Windows, use environment-win.yml instead of environment-dev.yml

Supported Python Versions
-------------------------

Python 3.6, 3.7 and 3.8 are officially supported, including testing during
development and packaging. Support for Python 2.7 has been dropped as of
August 6, 2019. Other Python versions, such as 3.9 and 3.5 and older, may
successfully build and function but no guarantee is made.

Testing your installation
-------------------------

mBuild uses `py.test <https://docs.pytest.org/en/stable/>`_ for unit testing. To run them simply run the following while in the base directory::

    $ conda install pytest
    $ py.test -v

Building the documentation
--------------------------

mBuild uses `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ to build its documentation. To build the docs locally, run the following while in the ``docs`` directory::
    
    $ pip install -r requirements.txt
    $ make html
