import glob
import os

from six import string_types

from mbuild.utils.io import load_top_opls, resource_filename
from mbuild.tools.parameterize.forcefield import prepare_atoms
from mbuild.tools.parameterize.atomtyper import find_atomtypes
from mbuild.tests.base_test import BaseTest


class TestOPLS(BaseTest):

    def test_all_molecules(self, only_run=None):
        resource_dir = resource_filename('mbuild', '../opls_validation')
        top_files = glob.glob(os.path.join(resource_dir, '*.top'))

        # Please update this file if you implement atom typing for a test case.
        implemented_tests_path = os.path.join(resource_dir, 'implemented_opls_tests.txt')
        correctly_implemented = [line.strip() for line in open(implemented_tests_path)]

        for top in top_files:
            top_name = os.path.split(top)[-1]
            loaded = load_top_opls(top)
            if loaded:
                topology, known_opls_types, mol_name = loaded
            else:
                continue
            if only_run and only_run != mol_name:
                continue
            elif mol_name not in correctly_implemented:
                continue
            print("Typing {} ({})...".format(mol_name, top_name))
            prepare_atoms(topology)
            find_atomtypes(topology._atoms, forcefield='OPLS-AA', debug=False)

            generated_opls_types = list()

            for i, atom in enumerate(topology._atoms):
                message = ('Found multiple or no OPLS types for atom {} in {} ({}): {}\n'
                           'Should be atomtype: {}'.format(
                    i, mol_name, top_name, atom.atomtype, known_opls_types[i]))
                assert isinstance(atom.atomtype, string_types), message
                generated_opls_types.append(atom.atomtype)
            both = zip(generated_opls_types, known_opls_types)
            message = "Found inconsistent OPLS types in {} ({}): {}".format(
                mol_name, top_name, list(zip(range(len(generated_opls_types)),
                                        generated_opls_types,
                                        known_opls_types)))

            assert all([a == b for a, b in both]), message
            print("Passed.\n")


if __name__ == "__main__":
    import pdb
    # TestOPLS().test_all_molecules('benzene')
    TestOPLS().test_all_molecules()
