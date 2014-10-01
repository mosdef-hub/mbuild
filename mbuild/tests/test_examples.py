#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.examples` module. """
import os

import mbuild.examples


class TestExamples:

    def test_all(self):
        sub_dirs = [os.path.basename(x[0]) for x in os.walk(os.path.dirname(mbuild.examples.__file__))][1:]

        for name in sub_dirs:
            example = __import__('mbuild.examples.{0}.{0}'.format(name), fromlist=[''])
            example.main()

