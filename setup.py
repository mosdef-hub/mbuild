"""mBuild: A hierarchical, component based molecule builder

With just a few lines of mBuild code, you can assemble reusable components into
complex molecular systems for molecular dynamics simulations. mBuild is
designed to minimize or even eliminate the need to explicitly translate and
orient components when building systems: you simply tell it to connect two
pieces! mBuild also keeps track of the system's topology so you don't have to
worry about manually defining bonds when constructing chemically bonded
structures from smaller components.
"""

import os
import subprocess
from pathlib import Path
from setuptools import setup, find_packages
from distutils.spawn import find_executable

#####################################
NAME = 'mbuild'
VERSION = "0.10.5"
ISRELEASED = True
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + '.dev0'
#####################################


def git_version():
    # Return the git revision as a string
    # copied from numpy setup.py
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = 'Unknown'

    return GIT_REVISION


def write_version_py(version, isreleased, filename):
    cnt = """
# This file is generated in setup.py at build time.
version = '{version}'
short_version = '{short_version}'
full_version = '{full_version}'
git_revision = '{git_revision}'
release = {release}
"""
    base_path = Path(__file__).parent
    file_path = (base_path / filename)
    # git_revision
    if os.path.exists('.git'):
        git_revision = git_version()
    else:
        git_revision = 'Unknown'

    # short_version, full_version
    if isreleased:
        full_version = version
        short_version = version
    else:
        full_version = ("{version}+{git_revision}"
                        .format(version=version, git_revision=git_revision))
        short_version = version

    with file_path.open('w') as f:
        f.write(cnt.format(version=version,
                           short_version=short_version,
                           full_version=full_version,
                           git_revision=git_revision,
                           release=isreleased))

def proto_procedure():
    # Find the Protocol Compiler and compile protocol buffers
    # Taken from https://github.com/protocolbuffers/protobuf/blob/fcfc47d405113b59bd43c2e54daf5d9fe5c44593/python/setup.py
    # Only compile if a protocompiler is found, otherwise don't do anything
    if 'PROTOC' in os.environ and os.path.exists(os.environ['PROTOC']):
      protoc = os.environ['PROTOC']
    elif os.path.exists("../src/protoc"):
      protoc = "../src/protoc"
    elif os.path.exists("../src/protoc.exe"):
      protoc = "../src/protoc.exe"
    elif os.path.exists("../vsprojects/Debug/protoc.exe"):
      protoc = "../vsprojects/Debug/protoc.exe"
    elif os.path.exists("../vsprojects/Release/protoc.exe"):
      protoc = "../vsprojects/Release/protoc.exe"
    else:
      protoc = find_executable("protoc")
      if protoc is None:
          protoc = find_executable("protoc.exe")

    if protoc is not None:
        compile_proto(protoc)


def compile_proto(protoc):
    protoc_command = [protoc, '-I=mbuild/formats/', 
            '--python_out=mbuild/formats/', 'compound.proto']
    subprocess.call(protoc_command)


if __name__ == '__main__':

    write_version_py(VERSION, ISRELEASED, 'mbuild/version.py')

    proto_procedure()

    setup(
        name=NAME,
        version=__version__,
        description=__doc__.split('\n'),
        long_description=__doc__,
        author='Janos Sallai, Christoph Klein',
        author_email='janos.sallai@vanderbilt.edu, christoph.klein@vanderbilt.edu',
        url='https://github.com/mosdef-hub/mbuild',
        download_url='https://github.com/mosdef-hub/mbuild/tarball/{}'.format(__version__),
        packages=find_packages(),
        package_data={'mbuild': ['utils/reference/*.{pdb,mol2}',
                                 'lib/*.{pdb,mol2}',
                                 ]},
        entry_points={
            'mbuild.plugins':[
                "Alkane = mbuild.lib.recipes.alkane:Alkane",
                "Monolayer = mbuild.lib.recipes.monolayer:Monolayer",
                "Polymer = mbuild.lib.recipes.polymer:Polymer",
                "SilicaInterface = mbuild.lib.recipes.silica_interface:SilicaInterface",
                "TiledCompound = mbuild.lib.recipes.tiled_compound:TiledCompound",
            ]
        },
        package_dir={'mbuild': 'mbuild'},
        include_package_data=True,
        license="MIT",
        zip_safe=False,
        keywords='mbuild',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Operating System :: MacOS',
        ],
    )

