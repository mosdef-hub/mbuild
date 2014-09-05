#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.hoomdxmlfile` module. """

from mbuild.examples.ethane.ethane import Ethane
from mbuild.testing.tools import get_fn


class TestHoomdXml:

    # def test_load_and_create(self):
    #     methyl = load_mol2('methyl.mol2')
    #
    def test_write(self):
        ethane = Ethane()
        ethane.save("ethane.hoomdxml")

    def test_update_from_file(self):
        ethane = Ethane()
        ethane.update_from_file(get_fn("ethane.hoomdxml"))

if __name__=="__main__":
    TestHoomdXml().test_write()
    TestHoomdXml().test_update_from_file()
