#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.formats.xyz` module. """
import pytest


class TestXyz:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def ethane(self):
        from mbuild.examples.ethane.ethane import Ethane
        return Ethane()

    def test_save(self, ethane):
        ethane.save(filename='ethane_out.xyz')


