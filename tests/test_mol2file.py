#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.mol2file` module. """

from mbuild.file_formats.mol2file import load_mol2, write_mol2
from mbuild.compound import Compound

class TestMol2:

    def test_load_and_create(self):
        methyl = load_mol2('methyl.mol2')

    def test_load_into(self):
        methyl = Compound()
        load_mol2('methyl.mol2', component=methyl)

    def test_write(self):
        methyl = load_mol2('methyl.mol2')
        write_mol2(methyl, filename='methyl_out.mol2')
