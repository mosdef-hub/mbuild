#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbuild.tools.parameterize` module. """

import glob
import os

from opls_validation.gromacs import load_top
from mbuild.tools.parameterize.oplsaa import opls_atomtypes


class TestTools:

    def test_all_molecules(self, only_run=None):
        top_files = glob.glob('../../opls_validation/*.top')
        for top in top_files[::-1]:
            top_name = os.path.split(top)[-1]
            loaded = load_top(top)
            if loaded:
                compound, known_opls_types, mol_name = loaded
            else:
                continue

            if only_run and only_run != mol_name:
                continue

            print "Typing {} ({})...".format(mol_name, top_name)

            opls_atomtypes(compound)
            generated_opls_types = list()
            for i, atom in enumerate(compound.atoms):
                # TODO: make less dirty
                opls_type = atom.opls_whitelist - atom.opls_blacklist
                opls_type = [a for a in opls_type]

                message = "Found multiple or no OPLS types for atom {} in {} ({}): {}".format(
                        i, mol_name, top_name, opls_type)
                assert len(opls_type) == 1, message
                generated_opls_types.append(opls_type[0])

            both = zip(generated_opls_types, known_opls_types)
            message = "Found inconsistent OPLS types in {} ({}): {}".format(
                    mol_name, top_name, zip(generated_opls_types, known_opls_types))
            assert all([a == b for a, b in both]), message
            print "passed.\n"


if __name__ == "__main__":
    import pdb
    TestTools().test_all_molecules()