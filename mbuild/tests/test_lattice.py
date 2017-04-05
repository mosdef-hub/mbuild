import numpy as np
import pytest
from mbuild.tests.base_test import BaseTest
import mbuild as mb


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
            neg_lattice = mb.Lattice(neg_test2, dimension=2,
                                     lattice_vectors=None)
        with pytest.raises(ValueError):
            zero_lattice = mb.Lattice(zero_test3, dimension=3,
                                      lattice_vectors=None)
        with pytest.raises(ValueError):
            neg_lattice = mb.Lattice(neg_test3, dimension=3,
                                     lattice_vectors=None)

    def test_lattice_spacings_correct(self):
        correct_spacings = [0.12, 0.13, 0.14]
        input_spacings = [0.12, 0.13, 0.14]
        test_lattice = mb.Lattice(input_spacings)
        np.testing.assert_array_equal(test_lattice.lattice_spacings,
                                      correct_spacings)

    def test_lattice_basis_atoms_default(self):
        space_1d = [1.]
        space_2d = [1., 1.]
        space_3d = [1., 1., 1.]
        correct_name = 'default'

        test_1d = mb.Lattice(space_1d, dimension=1, basis_atoms=None)
        test_2d = mb.Lattice(space_2d, dimension=2, basis_atoms=None)
        test_3d = mb.Lattice(space_3d, dimension=3, basis_atoms=None)

        assert len(test_1d.basis_atoms.keys()) == 1
        assert len(test_2d.basis_atoms.keys()) == 1
        assert len(test_3d.basis_atoms.keys()) == 1
        for lattice in test_1d, test_2d, test_3d:
            for key in lattice.basis_atoms.keys():
                assert key == correct_name

            correct_pos = [0. for x in range(lattice.dimension)]
            np.testing.assert_array_equal(lattice.basis_atoms[correct_name][0],
                                          correct_pos)

    def test_lattice_basis_atoms_input_type(self):
        incorrect_types = [list(), tuple(), str()]
        lat_space = [1., 1., 1.]
        for the_type in incorrect_types:
            with pytest.raises(TypeError):
                test_lat = mb.Lattice(lat_space, basis_atoms=the_type)

    def test_lattice_basis_atoms_input_size(self):
        dim = 3
        lat_space = [1., 1., 1.]
        basis_too_large = {'test': [[0., 0., 0., 0.]]}
        basis_too_large_mult = {'test': [[0., 0., 0., 0.], [0., 0., 0]]}
        basis_too_small = {'test': [[0., 0.]]}
        basis_too_small_mult = {'test': [[0., 0.], [0., 0., 0.]]}
        basis_empty = {'test': [[]]}
        basis_none = {'test': [[None, None, None]]}
        to_test = [basis_too_large, basis_too_small, basis_empty, basis_none]

        for errors in to_test:
            with pytest.raises(ValueError):
                test_lattice = mb.Lattice(lat_space, dimension=3,
                                          basis_atoms=errors)

    def test_lattice_basis_atoms_location(self):
        dim = 3
        lat_space = [1., 1., 1.]
        basis_neg = {'test': [[0., 0., -.2]]}
        basis_neg_mult = {'test': [[0., 0., 0.], [0., 0., -.2]]}
        basis_positive = {'test': [[0., 0., 1.2]]}
        basis_positive_mult = {'test': [[0., 0., 0.], [0., 0., 1.2]]}
        to_test = [basis_neg, basis_neg_mult, basis_positive, basis_positive_mult]

        for errors in to_test:
            with pytest.raises(ValueError):
                test_lattice = mb.Lattice(lat_space, dimension=dim, basis_atoms=errors)

    def test_lattice_basis_atoms_overlap(self):
        dim = 3
        lat_space = [1., 1., 1.]
        overlap_same_basis_repeat = {'test': [[0., 0., 0.], [0., 0., 0.]]}
        overlap_diff_basis_repeat = {'test': [[0., 0., 0.]],
                                     'test2': [[.5, .5, .5], [0., 0., 0.]]}

        to_test = [overlap_same_basis_repeat, overlap_diff_basis_repeat]
        for errors in to_test:
            with pytest.raises(ValueError):
                test_lattice = mb.Lattice(lat_space, dimension=dim,
                                          basis_atoms=errors)

    def test_lattice_overdefinied_inputs(self):
        dim = 3
        lat_space = [1., 1., 1.]
        angles = [90., 90., 90.]
        lattice_vec = [1, 1, 1]
