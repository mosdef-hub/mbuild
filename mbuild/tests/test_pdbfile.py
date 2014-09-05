#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.pdbfile` module. """

from mbuild.compound import Compound
from mbuild.examples.ethane.methyl import Methyl
from mbuild.testing.tools import get_fn


class TestPdb:

    def test_load_and_create(self):
        methyl = Compound.load(get_fn('methyl.pdb'))

    def test_update_from_file(self):
        methyl = Methyl()
        methyl.update_from_file(get_fn("methyl.pdb"))

    def test_write(self):
        methyl = Compound.load(get_fn('methyl.pdb'))
        methyl.save(filename='methyl_out.pdb')

if __name__=="__main__":
    TestPdb().test_load_and_create()
    TestPdb().test_update_from_file()
    TestPdb().test_write()

