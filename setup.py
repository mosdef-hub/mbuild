"""mBuild: A hierarchical, component based molecule builder

With just a few lines of mBuild code, you can assemble reusable components into
complex molecular systems for molecular dynamics simulations. mBuild is
designed to minimize or even eliminate the need to explicitly translate and
orient components when building systems: you simply tell it to connect two
pieces! mBuild also keeps track of the system's topology so you don't have to
worry about manually defining bonds when constructing chemically bonded
structures from smaller components.
"""

from __future__ import print_function

import os
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

import mbuild.version

try:
    import mdtraj
except ImportError:
    print('Building and running mbuild requires mdtraj. See '
          'http://mdtraj.org/latest/installation.html for help!', file=sys.stderr)
    sys.exit(1)

try:
    import scipy
except ImportError:
    print('Building and running mbuild requires scipy. We '
          'recommend using conda to install this: '
          'https://store.continuum.io/cshop/anaconda/', file=sys.stderr)
    sys.exit(1)

requirements = [line.strip() for line in open('requirements.txt').readlines()]

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(['mbuild'])
        sys.exit(errcode)

setup(
    name='mbuild',
    version=mbuild.version.short_version,
    description=__doc__.split('\n'),
    long_description=__doc__,
    author='Janos Sallai, Christoph Klein',
    author_email='janos.sallai@vanderbilt.edu, christoph.klein@vanderbilt.edu',
    url='https://github.com/sallai/mbuild',
    download_url='https://github.com/sallai/mbuild/tarball/{}'.format(
        mbuild.version.short_version),
    packages=find_packages(),
    package_data={'mbuild': ['utils/reference/*']},
    package_dir={'mbuild': 'mbuild'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT",
    zip_safe=False,
    keywords='mbuild',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License,'
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],
    test_suite='tests',
    cmdclass={'test': PyTest,
    },
    extras_require={'utils': ['pytest'],
    },
)
