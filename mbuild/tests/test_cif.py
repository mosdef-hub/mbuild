import numpy as np
import pytest
from operator import itemgetter
from collections import OrderedDict

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

    @pytest.mark.skipif(not has_garnett, reason="garnett package not installed")
    @pytest.mark.skipif(not has_pycifrw, reason="pycifrw package not installed")
    def test_cif_vs_manual(self):
        spacing = [0.760296570, 0.760296570, 0.437540800]
        points_dict = {"La": [[1/3, 2/3, 1/4],
                             [2/3, 1/3, 3/4]],
                       "Cl": [[0.69490400, 0.08690400, 1/4],
                              [0.60800000, 0.69490400, 3/4],
                              [0.30509600, 0.91309600, 3/4],
                              [0.39200000, 0.30509600, 1/4],
                              [0.91309600, 0.60800000, 1/4],
                              [0.08690400, 0.39200000, 3/4]]}

        lattice_manual = mb.Lattice(lattice_spacing=spacing, lattice_points=points_dict, angles=[90, 90, 120])
        lattice_cif = load_cif(file_or_path=get_fn("LaCl3.cif"))

        assert np.all(np.isclose(lattice_manual.lattice_spacing, lattice_cif.lattice_spacing))
        assert np.all(np.isclose(lattice_manual.angles, lattice_cif.angles))

        # sort dicts first (not necessary once we support py 3.7+ only)
        # dict sorted by keys
        dict_manual = OrderedDict(sorted(lattice_manual.lattice_points.items(), key=lambda t: t[0]))
        dict_cif = OrderedDict(sorted(lattice_cif.lattice_points.items(), key=lambda t: t[0]))
        keys_m = dict_manual.keys()
        keys_c = dict_cif.keys()

        for k_man, k_cif in zip(keys_m, keys_c):
            # sort the lists of lists
            points_man = dict_manual[k_man]
            points_cif = dict_cif[k_cif]
            points_man.sort()
            points_cif.sort()

            points_man = np.asarray(points_man)
            points_cif = np.asarray(points_cif)

            assert np.all(np.isclose(points_man, points_cif))
