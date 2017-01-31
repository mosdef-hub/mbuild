import numpy as np
import pytest
from mbuild.tests.base_test import BaseTest
import mbuild as mb
from lattice import Lattice
from collections import defaultdict


class TestLattice(BaseTest):
    """
    Unit Tests for Lattice class functionality.
    """
    def test_dimension_default(self):
        a_test = Lattice(dimension=None)
        assert a_test.dimension == 3

    def test_dimension_2D(self):
        a_test = Lattice(dimension=2)
        assert a_test.dimension == 2

    def test_dimension_3D(self):
        a_test = Lattice(dimension=3.0)
        assert a_test.dimension == 3

    def test_invalid_dimensions(self):
        with pytest.raises(ValueError):
            a_test = Lattice(dimension=1)
            a_test = Lattice(dimension=4)
        with pytest.raises(TypeError):
            a_test = Lattice(dimension='3')
            a_test = Lattice(dimension=([1, 2, 3]))

    def test_lattice_vectors_default(self):
        # default behavior for 2D and 3D
        two_dim_default = np.asarray(([1.0, 0.0], [0.0, 1.0]), dtype=float)
        three_dim_default = np.asarray(([1.0, 0.0, 0.0],
                                        [0.0, 1.0, 0.0],
                                        [0.0, 0.0, 1.0]), dtype=float)

        two_d_lattice = Lattice(dimension=2, lattice_vectors=None)
        three_d_lattice = Lattice(dimension=3, lattice_vectors=None)

        np.testing.assert_array_equal(two_dim_default,
                                      two_d_lattice.lattice_vectors)
        np.testing.assert_array_equal(three_dim_default,
                                      three_d_lattice.lattice_vectors)

    def test_lattice_vectors_invalid_shape(self):
        invalid_2d = np.asarray(([1, 0, 0], [0, 1, 0], [0, 0, 1]), dtype=float)
        invalid_3d = np.asarray(([1, 0], [0, 1]), dtype=float)
        with pytest.raises(ValueError):
            a_test_2d = Lattice(dimension=2, lattice_vectors=invalid_2d)
            a_test_3d = Lattice(dimension=3, lattice_vectors=invalid_3d)

    def test_colinear_lattice_vectors(self):
        invalid_2d = np.asarray(([1, 0], [3, 0]), dtype=float)
        invalid_3d = np.asarray(([1, 0, 0], [0, 1, 0], [2, 0, 0]), dtype=float)
        with pytest.raises(ValueError):
            a_test_2d = Lattice(dimension=2, lattice_vectors=invalid_2d)
            a_test_3d = Lattice(dimension=3, lattice_vectors=invalid_3d)

    def test_handedness_lattice_vectors(self):
        invalid_2d = np.asarray(([1, 2], [2, 1]), dtype=float)
        invalid_3d = np.asarray(([1, 2, 3], [3, 2, 1], [2, 1, 3]), dtype=float)
        with pytest.raises(ValueError):
            a_test_2d = Lattice(dimension=2, lattice_vectors=invalid_2d)
            a_test_3d = Lattice(dimension=3, lattice_vectors=invalid_3d)

    def test_lattice_spacings_default(self):
        spacing_test = Lattice(dimension=2, lattice_vectors=None,
                               lattice_spacings=None)
        np.testing.assert_array_equal(spacing_test.lattice_spacings,
                                      np.asarray([1, 1], dtype=float))
    def test_lattice_spacings_dimension(self):
        with pytest.raises(ValueError):
            spacing_test = Lattice(dimension=3, lattice_vectors=None,
                                   lattice_spacings=([.12], [.13], [.14]))
            spacing_test = Lattice(dimension=3, lattice_vectors=None,
                                   lattice_spacings=([.12, .13, .14, .15]))

    def test_lattice_spacings_negative_or_zero(self):
        with pytest.raises(ValueError):
            zero_test = ([.12, 0, .13])
            neg_test = ([.12, .13, -.14])
            zero_lattice = Lattice(dimension=3, lattice_vectors=None,
                                   lattice_spacings=zero_test)
            neg_lattice = Lattice(dimension=3, lattice_vectors=None,
                                  lattice_spacings=neg_test)

    def test_basis_default(self):
        vec_3D = defaultdict(list)
        vec_2D = defaultdict(list)
        vec_1D = defaultdict(list)
