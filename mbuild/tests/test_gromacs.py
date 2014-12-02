#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.formats.gromacs` module. """
import pytest


class TestGromacs:

    @pytest.fixture
    def ethane(self):
        from mbuild.examples.ethane.ethane import Ethane
        return Ethane()

    def test_save(self, ethane):
        ethane.save(filename='ethane_out.gro')


