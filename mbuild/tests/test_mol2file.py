#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.formats.mol2` module. """
from mbuild.examples.ethane.methyl import Methyl
from mbuild.trajectory import Trajectory
from mbuild.formats.mol2 import save_mol2
from mbuild.testing.tools import get_fn


class TestMol2:

    def test_load_and_create(self):
        methyl = Trajectory.load(get_fn('methyl.mol2'))
        methyl.to_compound()

    def test_load_into(self):
        methyl = Methyl()
        methyl.update_from_file(get_fn('methyl.mol2'))

    def test_save(self):
        methyl = Trajectory.load(get_fn('methyl.mol2'))
        save_mol2(methyl, filename='methyl_out.mol2')
