#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.pdbfile` module. """

from mbuild.file_formats.pdbfile import load_pdb, write_pdb
from mbuild.compound import Compound

class TestPdb:

    def test_load_and_create(self):
        methyl = load_pdb('methyl.pdb')

    def test_load_into(self):
        methyl = Compound()
        load_pdb('methyl.pdb', component=methyl)

    def test_write(self):
        methyl = load_pdb('methyl.pdb')
        write_pdb(methyl, filename='methyl_out.pdb') 
