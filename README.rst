mBuild: a hierarchical, component based molecule builder
========================================================

.. image:: https://badge.fury.io/py/mbuild.png
    :target: http://badge.fury.io/py/mbuild
    
.. image:: https://travis-ci.org/sallai/mbuild.png?branch=develop
        :target: https://travis-ci.org/sallai/mbuild
        
.. image:: https://coveralls.io/repos/sallai/mbuild/badge.png?branch=develop 
        :target: https://coveralls.io/r/sallai/mbuild?branch=develop

.. image:: https://readthedocs.org/projects/mbuild/badge/?version=develop
        :target: http://mbuild.readthedocs.org/en/develop/
        :alt: Documentation Status

mBuild is a component based molecule builder tool used to assemble complex
molecular systems from reusable parts for molecular dynamics simulations.

* Documentation: http://mbuild.rtfd.org/en/master/

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
        exit point for molecule data.

Optional packages:

    `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_
        VMD (Visual Molecular Dynamics) is a visualization program. Currently,
        we only use a (very) thin wrapper to call it from the command-line
        when you want to look at a Compound that you've built. A more robust
        integration with partial control from Python is planned in the near
        future (see issue [#32](/../../issues/32)).

To make your life easier, we recommend that you use a pre-packaged Python
distribution like `Enthought's Canopy <https://www.enthought.com/products/canopy/>`_
or `Continuum's Anaconda <https://store.continuum.io/>`_ in order to get all
of the dependencies.

=========================
Testing your installation
=========================

mBuild uses `py.test` for unit testing. To run them simply type run the
following while in the base directory:

 $ py.test

We need a LOT more tests so any help here is especially welcome!

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given. 

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/sallai/mbuild/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

mbuild could always use more documentation, whether as part of the 
official mbuild docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/sallai/mbuild/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `mbuild` for local development.

1. Fork the `mbuild` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/mbuild.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv mbuild
    $ cd mbuild/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature
   
   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests::

    $ flake8 mbuild tests
    $ python setup.py test

   To get flake8, just pip install it into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring.
3. The pull request should work for Python 2.7 (we intend to provide support for 2.6 and 3.3+ in the near future). Check 
   https://travis-ci.org/sallai/mbuild/pull_requests
   and make sure that the tests pass for all supported Python versions.

