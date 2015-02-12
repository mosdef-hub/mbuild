#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.formats.mol2` module. """
import pytest
from mbuild.components.small_groups.ch3 import Ch3
from mbuild.trajectory import Trajectory
from mbuild.testing.tools import get_fn

from base_test import BaseTest

class TestMol2(BaseTest):

    def test_load_and_create(self):
        methyl = Trajectory.load(get_fn('methyl.mol2'))
        methyl.to_compound()

    def test_load_into(self):
        methyl = Ch3()
        methyl.update_from_file(get_fn('methyl.mol2'))

    def test_save(self):
        methyl = Trajectory.load(get_fn('methyl.mol2'))
        methyl.save(filename='methyl_out.mol2')
