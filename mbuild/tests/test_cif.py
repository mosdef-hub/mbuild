import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.lattice import load_cif
from mbuild.utils.io import get_fn, has_garnett, has_pycifrw


class TestCif(BaseTest):
    """
    Unit tests for CIF file loading and Lattice generation.
    """
    @pytest.mark.skipif(not has_garnett, reason="garnett package not installed")
    @pytest.mark.skipif(not has_pycifrw, reason="pycifrw package not installed")
    def test_malformed_cif(self):
        with pytest.raises(Exception):
            load_cif(file_or_path=get_fn("extra_blank_field.cif"))

    @pytest.mark.skipif(not has_garnett, reason="garnett package not installed")
    @pytest.mark.skipif(not has_pycifrw, reason="pycifrw package not installed")
    def test_wrap_false(self):
        with pytest.raises(ValueError):
            load_cif(file_or_path=get_fn("needs_to_be_wrapped.cif"), wrap_coords=False)

    @pytest.mark.skipif(not has_garnett, reason="garnett package not installed")
    @pytest.mark.skipif(not has_pycifrw, reason="pycifrw package not installed")
    def test_wrap_true(self):
        assert load_cif(file_or_path=get_fn("needs_to_be_wrapped.cif"), wrap_coords=True)
