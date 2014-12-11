#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild` module. """
import pytest

class TestMbuild:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    def setUp(self):
        pass

    def test_something(self):
        pass

    def tearDown(self):
        pass
