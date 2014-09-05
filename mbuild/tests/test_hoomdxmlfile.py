#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.hoomdxmlfile` module. """

from mbuild.examples.ethane.ethane import Ethane
from mbuild.compound import Compound


class TestHoomdXml:

    # def test_load_and_create(self):
    #     methyl = load_mol2('methyl.mol2')
    #
    def test_load_into(self):
        ethane = Compound()
        ethane.update_from_file("ethane.hoomdxml", relative_to_module=__name__)

    def test_write(self):
        ethane = Ethane()
        ethane.save("ethane.hoomdxml")


if __name__=="__main__":
    TestHoomdXml().test_load_into()