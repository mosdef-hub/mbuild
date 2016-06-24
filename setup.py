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
doc_lines = __doc__.split('\n')

import os
import sys
from setuptools import setup, Extension, find_packages
sys.path.insert(0, '.')
from basesetup import (write_version_py, build_ext,
                       StaticLibrary, CompilerDetection)

try:
    import Cython
    if Cython.__version__ < '0.19':
        raise ImportError
    from Cython.Build import cythonize
except ImportError:
    print('-'*80, file=sys.stderr)
    print('''Error: building mbuild requires cython>=0.19
Try running the command ``pip install cython`` or
``conda install cython`` or see http://cython.org/ for more information.
If you're feeling lost, we recommend downloading the (free) Anaconda python
distribution https://www.continuum.io/downloads, because it comes with
these components included.''', file=sys.stderr)
    print('-'*80, file=sys.stderr)
    sys.exit(1)


try:
    # add an optional --disable-openmp to disable OpenMP support
    sys.argv.remove('--disable-openmp')
    disable_openmp = True
except ValueError:
    disable_openmp = False


#####################################
VERSION = "0.6.0.dev0"
ISRELEASED = False
__version__ = VERSION
write_version_py(VERSION, ISRELEASED, 'mbuild/version.py')
#####################################


# Global info about compiler
compiler = CompilerDetection(disable_openmp)
compiler.initialize()

extra_cpp_libraries = []
if sys.platform == 'darwin':
    extra_cpp_libraries.append('stdc++')
    os.environ['CXX'] = 'clang++'
    os.environ['CC'] = 'clang'
if sys.platform == 'win32':
    extra_cpp_libraries.append('Ws2_32')
    # For determining if a path is relative (for dtr)
    extra_cpp_libraries.append('Shlwapi')

setup(
    name='mbuild',
    author='Janos Sallai, Christoph Klein',
    author_email='janos.sallai@vanderbilt.edu, christoph.klein@vanderbilt.edu',
    description=doc_lines[0],
    long_description='\n'.join(doc_lines[2:]),
    version=__version__,
    license="MIT",
    url='https://github.com/imodels/mbuild',
    download_url='https://github.com/imodels/mbuild/tarball/{}'.format(__version__),
    platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
    keywords='mbuild',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],
    packages=find_packages(),
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize('mbuild/bond_graph.pyx'),
    package_data={'mbuild': ['utils/reference/*.{pdb,mol2}',
                             'lib/*.{pdb,mol2}',
                             ]},
    package_dir={'mbuild': 'mbuild'},
    include_package_data=True,
    zip_safe=False,
)
