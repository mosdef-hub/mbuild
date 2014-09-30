#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.compound` module. """

from mbuild.compound import Compound
from mbuild.examples.ethane.methyl import Methyl
from mbuild.testing.tools import get_fn


class TestCompound:

    def test_load_and_create(self):
        methyl = Compound.load(get_fn('methyl.pdb'))

    def test_update_from_file(self):
        methyl = Methyl()
        methyl.update_from_file(get_fn("methyl.pdb"))

    def test_write(self):
        methyl = Compound.load(get_fn('methyl.pdb'))
        methyl.save(filename='methyl_out.pdb')


