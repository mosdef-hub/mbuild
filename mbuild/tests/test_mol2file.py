#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.formats.mol2` module. """
import pytest
from mbuild.examples.ethane.methyl import Methyl
from mbuild.trajectory import Trajectory
from mbuild.testing.tools import get_fn


class TestMol2:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    def test_load_and_create(self):
        methyl = Trajectory.load(get_fn('methyl.mol2'))
        methyl.to_compound()

    def test_load_into(self):
        methyl = Methyl()
        methyl.update_from_file(get_fn('methyl.mol2'))

    def test_save(self):
        methyl = Trajectory.load(get_fn('methyl.mol2'))
        methyl.save(filename='methyl_out.mol2')
