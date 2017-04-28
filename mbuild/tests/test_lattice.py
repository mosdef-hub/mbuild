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

    def test_dimension_1d(self):
        space = [1, ]
        a_test = mb.Lattice(space, dimension=1)
        assert a_test.dimension == 1

    def test_dimension_2d(self):
        space = [1, 1]
        a_test = mb.Lattice(space, dimension=2)
        assert a_test.dimension == 2

    def test_dimension_3d(self):
        space = [1, 1, 1]
        a_test = mb.Lattice(space, dimension=3.0)
        assert a_test.dimension == 3

    @pytest.mark.parametrize("dim, space",
                            [
                                (4, [1, 1, 1, 1]),
                                (([1, 2, 3]), [1, 1, 1])
                            ]
                            )
    def test_invalid_dimensions(self, dim, space):
        with pytest.raises(Exception):
            test_lattice = mb.Lattice(space, dimension=dim)

    @pytest.mark.parametrize("spacing, vec_default, dimension",
                            [
                                ([1, ], np.identity(1, dtype=float), 1),
                                ([1, 1], np.identity(2, dtype=float), 2),
                                ([1, 1, 1], np.identity(3, dtype=float), 3)
                            ]
                            )
    def test_lattice_vectors_default(self, spacing, vec_default, dimension):

        test_lattice = mb.Lattice(spacing, dimension=dimension,
                                  lattice_vectors=None )
        np.testing.assert_array_equal(vec_default,
                                      test_lattice.lattice_vectors)

    @pytest.mark.parametrize("spacing, invalid_shape, dim",
                            [
                                ([1, ], np.identity(3, dtype=float), 1),
                                ([1, 1], np.identity(3, dtype=float), 2),
                                ([1, 1, 1], np.identity(2, dtype=float), 3)
                            ]
                            )
    def test_lattice_vectors_invalid_shape(self, spacing, invalid_shape, dim):
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(spacing, dimension=dim,
                                      lattice_vectors=invalid_shape)

    @pytest.mark.parametrize("space, dim, colinear",
                            [
                                ([1, 1], 2, [[1, 0], [3, 0]]),
                                ([1, 1, 1], 3, [[1, 0, 0], [0, 1, 0], [2, 0, 0]])
                            ]
                            )
    def test_colinear_lattice_vectors(self, space, dim, colinear):
        vectors = np.asarray(colinear, dtype=None)
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim,
                                      lattice_vectors=vectors)

    @pytest.mark.parametrize("spacing, dim, invalid",
                             [
                                ([1, 1], 2, [[1, 2], [2, 1]]),
                                ([1, 1, 1], 3, [[1, 2, 3], [3, 2, 1], [2, 1, 3]])
                             ]
                             )
    def test_handedness_lattice_vectors(self, spacing, dim, invalid):
        invalid_handed = np.asarray(invalid, dtype=float)
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(spacing, dimension=dim,
                                      lattice_vectors=invalid_handed)

    @pytest.mark.parametrize("dim, invalid_space",
                             [
                                (3, ([.12], [.13], [.14])),
                                (3, ([0.12, 0.13, 0.14, 0.15]))
                             ]
                             )
    def test_lattice_spacings_dimension(self, dim, invalid_space):
        with pytest.raises(ValueError):
            spacing_test = mb.Lattice(invalid_space, dimension=dim)

    @pytest.mark.parametrize("dim , space_invalid",
                             [
                                (1, [0]),
                                (1, [-0.14]),
                                (2, [.12, 0]),
                                (2, [.13, -.14]),
                                (3, [.12, 0, 0.13]),
                                (3, [.12, .13, -.14])
                             ]
                             )
    def test_lattice_spacings_negative_or_zero(self, dim , space_invalid):
        with pytest.raises(ValueError):
            zero_lattice = mb.Lattice(space_invalid, dimension=dim)

    def test_lattice_spacings_correct(self):
        correct_spacings = [0.12, 0.13, 0.14]
        input_spacings = [0.12, 0.13, 0.14]
        test_lattice = mb.Lattice(input_spacings)
        np.testing.assert_array_equal(test_lattice.lattice_spacings,
                                      correct_spacings)

    @pytest.mark.parametrize("space, correct, dim",
                            [
                                ([1.], 'default', 1),
                                ([1., 1.], 'default', 2),
                                ([1., 1., 1.], 'default', 3)
                            ]
                            )
    def test_basis_atoms_default(self, space, correct, dim):
        test_lattice = mb.Lattice(space, dimension=dim)
        assert len(test_lattice.basis_atoms.keys()) == 1

        for key in test_lattice.basis_atoms.keys():
            assert key == correct

        correct_pos = [0. for x in range(dim)]
        np.testing.assert_array_equal(test_lattice.basis_atoms[correct][0], correct_pos)

    @pytest.mark.parametrize("space, dim, type",
                             [
                                 ([1, 1, 1], 3, list()),
                                 ([1, 1, 1], 3, tuple()),
                                 ([1, 1, 1], 3, str())
                             ]
                             )
    def test_basis_atoms_input_type(self, space, dim, type):
        with pytest.raises(TypeError):
            test_lat = mb.Lattice(space, dimension=dim, basis_atoms=type)

    @pytest.mark.parametrize("space, dim, name, coords",
                             [
                                 ([1, 1, 1], 3, 'test', [[0, 0, 0, 0]]),
                                 ([1, 1, 1], 3, 'test', [[0, 0, 0], [0, 0, 0, 0]]),
                                 ([1, 1, 1], 3, 'test', [[0, 0], [0, 0, 0]]),
                                 ([1, 1, 1], 3, 'test', [[0, 0]])
                             ]
                             )
    def test_basis_atoms_input_size(self, space, dim, name, coords):
        basis = {name: coords}
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim, basis_atoms=basis)

    @pytest.mark.parametrize("space, dim, name, coords",
                             [
                                 ([1, 1, 1], 3, 'test', [[0., 0., -.2]]),
                                 ([1, 1, 1], 3, 'test', [[0., 0., 0.], [0., 0., -.2]]),
                                 ([1, 1, 1], 3, 'test', [[0., 0., 1.2]]),
                                 ([1, 1, 1], 3, 'test', [[0., 0., 0.], [0., 0., 1.2]])
                             ]
                             )
    def test_basis_atoms_location(self, space, dim, name, coords):
        basis = {name: coords}
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim, basis_atoms=basis)

    @pytest.mark.parametrize("space, dim, name, overlap",
                             [
                                 ([1, 1, 1], 3, 'test', [[0., 0., 0.], [0., 0., 0.]]),
                                 ([1, 1, 1], 3, 'test', [[0., .5, 0.], [0., .5, 0]])
                             ]
                             )
    def test_basis_overlap_1_tag(self, space, dim, name, overlap):
        basis = {name: overlap}
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim, basis_atoms=basis)

    def test_basis_atoms_overlap_2tags(self):
        dim = 3
        lat_space = [1., 1., 1.]
        overlap_diff_basis_repeat = {'test': [[0., 0., 0.]],
                                     'test2': [[.5, .5, .5], [0., 0., 0.]]}

        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(lat_space, dimension=dim,
                                      basis_atoms=overlap_diff_basis_repeat)

    def test_overdefinied_inputs(self):
        dim = 3
        lat_space = [1., 1., 1.]
        angles = [90., 90., 90.]
        lattice_vec = [1, 1, 1]

        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(lat_space, dimension=dim, angles=angles,
                                      lattice_vectors=lattice_vec)

    @pytest.mark.parametrize("dim, space, angles",
                             [
                                 (1, [1], [90]),
                                 (2, [1, 1], [90, 90]),
                                 (3, [1, 1, 1], [90, 90]),
                                 (2, [1, 1], [90, 90, 90]),
                                 (3, [1, 1, 1], [90, 90, 90, 90])
                             ]
                             )
    def test_angle_invalid_size(self, dim, space, angles):
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim, angles=angles)

    @pytest.mark.parametrize("dim, space, angle",
                             [
                                 (2, [1, 1], [180]),
                                 (2, [1, 1], [0]),
                                 (3, [1, 1, 1], [93, 70, 180]),
                                 (3, [1, 1, 1], [93, 70, 0])
                             ]
                             )
    def test_angle_too_large_small(self, dim, space, angle):
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim, angles=angle)

    @pytest.mark.parametrize("dim, space, angle",
                             [
                                 (3, [1, 1, 1], [90, 90, 185]),
                                 (3, [1, 1, 1], [-90, -90, -185])
                             ]
                             )
    def test_3d_sum_incorrect(self, dim, space, angle):
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim, angles=angle)

    @pytest.mark.parametrize("dim, space, angle",
                             [
                                 (3, [1, 1, 1], [60, 42, 105]),
                                 (3, [1, 1, 1], [65, 120, 42])
                             ]
                             )
    def test_each_angle_sum_incorrect(self, dim, space, angle):
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim, angles=angle)

    @pytest.mark.parametrize("angle",
                             [
                                 ([0]),
                                 ([180]),
                                 ([181]),
                                 ([-181])
                             ]
                             )
    def test_angle_2d_incorrect(self, angle):
        dim = 2
        space = [1, 1]
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(space, dimension=dim, angles=angle)

    @pytest.mark.parametrize("spacing, dim, x, y, z",
                               [
                                    ([1, 1, 1], 3, 0, 1, 1),
                                    ([1, 1, 1], 3, 1, 0, 1),
                                    ([1, 1, 1], 3, 1, 1, 0),
                                    ([1, 1, 1], 3, None, 1, 0),
                                    ([1, 1, 1], 3, -1, -1, -1)
                               ]
                             )
    def test_populate_3d_incorrect_inputs(self, spacing, dim, x, y, z):
        test_lattice = mb.Lattice(spacing, dimension=dim)
        with pytest.raises(ValueError):
            test_lattice.populate(x=x, y=y, z=z)

    @pytest.mark.parametrize("spacing, dim, x, y",
                             [
                                    ([1, 1], 2, 0, 1),
                                    ([1, 1], 2, 1, 0),
                                    ([1, 1], 2, None, 0),
                                    ([1, 1], 2, -1, -1)
                             ]
                             )
    def test_populate_2d_incorrect_inputs(self, spacing, dim, x, y):
        test_lattice = mb.Lattice(spacing, dimension=dim)
        with pytest.raises(ValueError):
            test_lattice.populate(x=x, y=y)

    @pytest.mark.parametrize("spacing, dim, x",
                             [
                                    ([1], 1, 0),
                                    ([1], 1, -1)
                             ]
                             )
    def test_populate_1d_incorrect_inputs(self, spacing, dim, x):
        test_lattice = mb.Lattice(spacing, dimension=dim)
        with pytest.raises(ValueError):
            test_lattice.populate(x=x)
            test_lattice.populate(y=x)
            test_lattice.populate(z=x)

    def test_populate_3d_default(self):
        test_lattice = mb.Lattice([1, 1, 1], dimension=3)
        a = test_lattice.populate()
        np.testing.assert_array_equal(a.xyz[0], [0., 0., 0.])

    def test_populate_2d_default(self):
        test_lattice = mb.Lattice([1, 1], dimension=2)
        a = test_lattice.populate()
        np.testing.assert_array_equal(a.xyz[0], [0., 0., 0.])

    def test_populate_1d_default(self):
        test_lattice = mb.Lattice([1], dimension=1)
        a = test_lattice.populate()
        np.testing.assert_array_equal(a.xyz[0], [0., 0., 0.])

    @pytest.mark.parametrize("my_type",
                             [
                                    ([]),
                                    (()),
                                    (np.array),
                                    (np.ndarray)
                             ]
                             )
    def test_populate_basis_type_incorrect(self, my_type):
        test_lattice = mb.Lattice([1, 1, 1], dimension=3)
        with pytest.raises(TypeError):
            test_lattice.populate(compound_dict=my_type)

    @pytest.mark.parametrize("not_compound",
                             [
                                    (1),
                                    (mb.Box(lengths=[1, 1, 1])),
                                    ("aLattice")
                             ]
                             )
    def test_populate_not_compound(self, not_compound):
        test_lattice = mb.Lattice([1, 1, 1])
        particle_dict = {'default': not_compound}
        with pytest.raises(TypeError):
            test_lattice.populate(compound_dict=particle_dict)
