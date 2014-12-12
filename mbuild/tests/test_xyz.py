#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.formats.xyz` module. """
import pytest
from base_test import BaseTest


class TestXyz(BaseTest):

    @pytest.fixture
    def ethane(self):
        from mbuild.examples.ethane.ethane import Ethane
        return Ethane()

    def test_save(self, ethane):
        ethane.save(filename='ethane_out.xyz')


