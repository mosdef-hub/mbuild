import numpy as np
import pytest
from mbuild.tests.base_test import BaseTest
import mbuild as mb
from collections import defaultdict


class TestLattice(BaseTest):
    """
    Unit Tests for Lattice class functionality.
    """
    def test_dimension_default(self):
        space = [1, 1, 1]
        a_test = mb.Lattice(dimension=None, lattice_spacings=space)
        assert a_test.dimension == 3
        a_test = mb.Lattice(space, dimension='3')
        assert a_test.dimension == 3

    def test_dimension_1D(self):
        space = [1, ]
        a_test = mb.Lattice(space, dimension=1)
        assert a_test.dimension == 1

    def test_dimension_2D(self):
        space = [1, 1]
        a_test = mb.Lattice(space, dimension=2)
        assert a_test.dimension == 2

    def test_dimension_3D(self):
        space = [1, 1, 1]
        a_test = mb.Lattice(space, dimension=3.0)
        assert a_test.dimension == 3

    def test_invalid_dimensions(self):
        space = [1, 1, 1, 1]
        space3 = [1, 1, 1]
        with pytest.raises(ValueError):
            a_test = mb.Lattice(space, dimension=4)
        with pytest.raises(TypeError):
            a_test = mb.Lattice(space3, dimension=([1, 2, 3]))

    def test_lattice_vectors_default(self):
        # default behavior for 2D and 3D
        space1 = [1, ]
        space2 = [1, 1]
        space3 = [1, 1, 1]
        one_dim_default = np.asarray(([1.0]), dtype=float)
        two_dim_default = np.asarray(([1.0, 0.0], [0.0, 1.0]), dtype=float)
        three_dim_default = np.asarray(([1.0, 0.0, 0.0],
                                        [0.0, 1.0, 0.0],
                                        [0.0, 0.0, 1.0]), dtype=float)
        one_d_lattice = mb.Lattice(space1, dimension=1, lattice_vectors=None)
        two_d_lattice = mb.Lattice(space2, dimension=2, lattice_vectors=None)
        three_d_lattice = mb.Lattice(space3, dimension=3, lattice_vectors=None)

        np.testing.assert_array_equal(one_dim_default,
                                      one_d_lattice.lattice_vectors)
        np.testing.assert_array_equal(two_dim_default,
                                      two_d_lattice.lattice_vectors)
        np.testing.assert_array_equal(three_dim_default,
                                      three_d_lattice.lattice_vectors)

    def test_lattice_vectors_invalid_shape(self):
        space1 = [1, ]
        space2 = [1, 1, ]
        space3 = [1, 1, 1]
        invalid_1d = np.asarray(([1, 0], [0, 1]), dtype=float)
        invalid_2d = np.asarray(([1, 0, 0], [0, 1, 0], [0, 0, 1]), dtype=float)
        invalid_3d = np.asarray(([1, 0], [0, 1]), dtype=float)
        with pytest.raises(ValueError):
            a_test_1d = mb.Lattice(space1, dimension=1,
                                lattice_vectors=invalid_1d)
        with pytest.raises(ValueError):
            a_test_2d = mb.Lattice(space2, dimension=2,
                                lattice_vectors=invalid_2d)
        with pytest.raises(ValueError):
            a_test_3d = mb.Lattice(space3, dimension=3,
                                lattice_vectors=invalid_3d)

    def test_colinear_lattice_vectors(self):
        shape2 = [1, 1]
        shape3 = [1, 1, 1]
        invalid_2d = np.asarray(([1, 0], [3, 0]), dtype=float)
        invalid_3d = np.asarray(([1, 0, 0], [0, 1, 0], [2, 0, 0]), dtype=float)
        with pytest.raises(ValueError):
            a_test_2d = mb.Lattice(shape2, dimension=2,
                                lattice_vectors=invalid_2d)
        with pytest.raises(ValueError):
            a_test_3d = mb.Lattice(shape3, dimension=3,
                                lattice_vectors=invalid_3d)

    def test_handedness_lattice_vectors(self):
        shape2 = [1, 1]
        shape3 = [1, 1, 1]
        invalid_2d = np.asarray(([1, 2], [2, 1]), dtype=float)
        invalid_3d = np.asarray(([1, 2, 3], [3, 2, 1], [2, 1, 3]), dtype=float)
        with pytest.raises(ValueError):
            a_test_2d = mb.Lattice(shape2, dimension=2,
                                lattice_vectors=invalid_2d)
        with pytest.raises(ValueError):
            a_test_3d = mb.Lattice(shape3, dimension=3,
                                lattice_vectors=invalid_3d)

    def test_lattice_spacings_dimension(self):
        with pytest.raises(ValueError):
            spacing_test = mb.Lattice(dimension=3, lattice_vectors=None,
                                   lattice_spacings=([.12], [.13], [.14]))

        with pytest.raises(ValueError):
            spacing_test = mb.Lattice(dimension=3, lattice_vectors=None,
                                   lattice_spacings=([.12, .13, .14, .15]))

    def test_lattice_spacings_negative_or_zero(self):
        zero_test1 = [0]
        neg_test1 = [-.14]
        zero_test2 = [.12, 0]
        neg_test2 = [.13, -.14]
        zero_test3 = [.12, 0, .13]
        neg_test3 = [.12, .13, -.14]
        with pytest.raises(ValueError):
            zero_lattice = mb.Lattice(zero_test1, dimension=1,
                                   lattice_vectors=None)
        with pytest.raises(ValueError):
            neg_lattice = mb.Lattice(neg_test1, dimension=1,
                                  lattice_vectors=None)
        with pytest.raises(ValueError):
            zero_lattice = mb.Lattice(zero_test2, dimension=2,
                                   lattice_vectors=None)
        with pytest.raises(ValueError):
            neg_lattice = mb.Lattice(neg_test2, dimension=2, lattice_vectors=None)
        with pytest.raises(ValueError):
            zero_lattice = mb.Lattice(zero_test3, dimension=3,
                                   lattice_vectors=None)
        with pytest.raises(ValueError):
            neg_lattice = mb.Lattice(neg_test3, dimension=3, lattice_vectors=None)

    def test_basis_default(self):
        three_d = mb.Lattice([1, 1, 1], dimension=3, basis_vectors=None)
        two_d = mb.Lattice([1, 1], dimension=2, basis_vectors=None)
        one_d = mb.Lattice([1], dimension=1, basis_vectors=None)

        assert len(three_d.basis_vectors) == 1
        assert len(two_d.basis_vectors) == 1
        assert len(one_d.basis_vectors) == 1

        assert three_d.basis_vectors['default'][0] == (0, 0, 0)
        assert two_d.basis_vectors['default'][0] == (0, 0)
        assert one_d.basis_vectors['default'][0] == (0, )

    def test_basis_1d(self):
        basis_1d = (('test', [.25]),)
        lat_1d = mb.Lattice([1], dimension=1, basis_vectors=basis_1d)
        assert len(lat_1d.basis_vectors) == 1
        assert lat_1d.basis_vectors['test'][0] == [.25]

    def test_basis_2d(self):
        basis_2d = (('test', [.25, .60]),)
        lat_1d = mb.Lattice([1, 1], dimension=2, basis_vectors=basis_2d)
        assert len(lat_1d.basis_vectors) == 1
        assert lat_1d.basis_vectors['test'][0] == [.25, .60]

    def test_basis_3d(self):
        basis_3d = (('test', [.25, .60, .70]),)
        lat_1d = mb.Lattice([1, 1, 1], dimension=3, basis_vectors=basis_3d)
        assert len(lat_1d.basis_vectors) == 1
        assert lat_1d.basis_vectors['test'][0] == [.25, .60, .70]

    def test_basis_multi(self):
        basis_1d = (('test1', [0]), ('test2', [.25]),)
        basis_2d = (('test1', [0, 0]), ('test2', [.25, .25]),)
        basis_3d = (('test1', [0, 0, 0]), ('test2', [.25, .25, .25]),)

        lat_1d = mb.Lattice([1], dimension=1, basis_vectors=basis_1d)
        lat_2d = mb.Lattice([1, 1], dimension=2, basis_vectors=basis_2d)
        lat_3d = mb.Lattice([1, 1, 1], dimension=3, basis_vectors=basis_3d)

        assert len(lat_1d.basis_vectors) == 2
        assert len(lat_2d.basis_vectors) == 2
        assert len(lat_3d.basis_vectors) == 2

        assert lat_1d.basis_vectors['test1'][0] == [0]
        assert lat_1d.basis_vectors['test2'][0] == [.25]
        assert lat_2d.basis_vectors['test1'][0] == [0, 0]
        assert lat_2d.basis_vectors['test2'][0] == [.25, .25]
        assert lat_3d.basis_vectors['test1'][0] == [0, 0, 0]
        assert lat_3d.basis_vectors['test2'][0] == [.25, .25, .25]
