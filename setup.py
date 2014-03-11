#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='mbuild',
    version='0.3.0',
    description='mbuild is a component based molecule builder tool to assemble large molecular systems from reusable parts for molecular dynamics simulations',
    long_description=readme + '\n\n' + history,
    author='Janos Sallai',
    author_email='janos.sallai@vanderbilt.edu',
    url='https://github.com/sallai/mbuild',
    packages=[
        'mbuild',
    ],
    package_dir={'mbuild': 'mbuild'},
    include_package_data=True,
    install_requires=[
    ],
    license="BSD",
    zip_safe=False,
    keywords='mbuild',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
)