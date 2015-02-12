#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.tools.parameterize` module. """

import glob
import os

from mbuild.testing import tools

from mbuild.testing.tools import load_top_opls
from mbuild.tools.parameterize.atomtyper import find_atomtypes
from base_test import BaseTest


class TestTools(BaseTest):
    def test_all_molecules(self, only_run=None):
        # Path doesn't work for py.test and is probably bad practice anyways.
        resource_dir = tools.resource_filename('mbuild', '../opls_validation')
        top_files = glob.glob(os.path.join(resource_dir, '*.top'))

        implemented_tests_path = os.path.join(resource_dir,
                                              'implemented_opls_tests.txt')
        correctly_implemented = [line.strip() for line in
                                 open(implemented_tests_path)]
        for top in top_files:
            top_name = os.path.split(top)[-1]
            loaded = load_top_opls(top)
            if loaded:
                compound, known_opls_types, mol_name = loaded
            else:
                continue
            if only_run and only_run != mol_name:
                continue
            elif mol_name not in correctly_implemented:
                continue


            print "Typing {} ({})...".format(mol_name, top_name)
            find_atomtypes(compound, forcefield='OPLS-AA', debug=True)
            generated_opls_types = list()
            for i, atom in enumerate(compound.atoms):
                message = ('Found multiple or no OPLS types for atom {} in {} ({}): {}\n'
                           'Should be atomtype: {}'.format(
                    i, mol_name, top_name, atom.atomtype, known_opls_types[i]))
                assert len(atom.atomtype) == 1, message
                generated_opls_types.append(atom.atomtype[0])

            both = zip(generated_opls_types, known_opls_types)
            message = "Found inconsistent OPLS types in {} ({}): {}".format(
                mol_name, top_name, zip(range(len(generated_opls_types)),
                                        generated_opls_types,
                                        known_opls_types))
            assert all([a == b for a, b in both]), message
            print "passed.\n"


if __name__ == "__main__":
    import pdb
    # TestTools().test_all_molecules('ethylene-carbonate')
    TestTools().test_all_molecules()
