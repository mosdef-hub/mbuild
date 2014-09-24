#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

from pip.req import parse_requirements

requirements_lines = [line.strip() for line in open('requirements.txt').readlines()]
reqs = list(filter(None, requirements_lines))

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
#        errcode = pytest.main(self.test_args)
        errcode = pytest.main(['mbuild'])
        sys.exit(errcode)

setup(
    name='mbuild',
    version='0.3.0',
    description='mbuild is a component based molecule builder tool to assemble large molecular systems from reusable parts for molecular dynamics simulations',
    long_description=readme + '\n\n' + history,
    author='Janos Sallai, Christoph Klein',
    author_email='janos.sallai@vanderbilt.edu, christoph.klein@vanderbilt.edu',
    url='https://github.com/sallai/mbuild',
    packages=[
        'mbuild',
    ],
    package_dir={'mbuild': 'mbuild'},
    include_package_data=True,
    install_requires=reqs,
    license="LGPL",
    zip_safe=False,
    keywords='mbuild',
    classifiers=[
        'Development Status :: 3 - Beta',
        'Intended Audience :: Computational Chemistry',
        'License :: OSI Approved :: LGPL License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    cmdclass={
      'test':PyTest,
    },
    extras_require={
      'testing': ['pytest'],
    },
)