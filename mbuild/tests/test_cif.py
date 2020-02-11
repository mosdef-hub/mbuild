import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.lattice import load_cif
from CifFile import StarError
from mbuild.utils.io import get_fn

#TODO: use the has_pycifrw importer

class TestCif(BaseTest):
    """
    Unit tests for CIF file loading and Lattice generation.
    """

    def test_malformed_cif(self):
        with pytest.raises(StarError):
            load_cif(file_or_path=get_fn("extra_blank_field.cif"))

    def test_has_spacings(self):
        pass

    def test_has_atomsites(self):
        #TODO Make sure that a cif file has atomsites
        pass
    
    def test_has_xyz_pos(self):
        #TODO Check that a cif file has positions for each atomsite
        pass
    
    def test_occupation_probability(self):
        #TODO Make sure the user is warned if a cif file has atomsites with probability distributions
        pass

    def test_infer_compounds(self):
        #TODO attempt to make compounds from atom_site_label
        pass
