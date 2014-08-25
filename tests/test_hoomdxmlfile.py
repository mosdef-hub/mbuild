#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.hoomdxmlfile` module. """

from mbuild.file_formats.hoomdxmlfile import write_hoomdxml
from mbuild.file_formats.mol2file import load_mol2

class TestMol2:

    def test_load_and_create(self):
        #methyl = load_hoomdxml('methyl.hoomdxml')
        pass

    def test_load_into(self):
        #methyl = Compound()
        #load_hoomdxml('methyl.hoomdxml', component=methyl)
        pass

    def test_write(self):
        methyl = load_mol2('methyl.mol2').
        write_hoomdxml(methyl, filename='methyl_out.xml')
