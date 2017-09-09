import numpy as np
import pytest
from mbuild.tests.base_test import BaseTest
import mbuild as mb


class TestLattice(BaseTest):
    """
    Unit Tests for Lattice class functionality.
    """

    @pytest.mark.parametrize("spacing",
                             [
                                ([1, 1, 1]),
                                ([0.1, 0.1, 0.1]),
                                (['1', '1', '1']),
                                (['1', 0.1, '0.1'])
                             ]
                             )
    def test_spacing_success(self, spacing):
        spacing = np.asarray(spacing, dtype=np.float64)
        spacing = np.reshape(spacing, (1, 3), order='C')
        test_lattice = mb.Lattice(lattice_spacing=spacing)
        np.testing.assert_allclose(spacing, test_lattice.lattice_spacing,
                                   rtol=1e-7, atol=0, equal_nan=True)

    @pytest.mark.parametrize("dim, spacing",
                             [
                                (3, [1, 1, 1]),
                                (3, [1, 1, 0]),
                                (3, [1, 0, 0])
                             ])
    def test_dimension_set(self, dim, spacing):
        test_lattice = mb.Lattice(lattice_spacing=spacing)
        assert test_lattice.dimension == dim

    @pytest.mark.parametrize("spacing",
                             [
                                ([1]),
                                (1),
                                ([1, 1]),
                                ([-1, 1, 1]),
                                ([1, 'a']),
                                ([]),
                                ([None, None, None]),
                             ])
    def test_spacing_incorrect(self, spacing):
        with pytest.raises(ValueError):
            mb.Lattice(lattice_spacing=spacing)

    @pytest.mark.parametrize("spacing",
                             [
                                ([0.1, 0.1, 0.1]),
                                ([1, 2, 3]),
                                (['1', '2', '3']),
                                ([1, 2, '3']),
                                ([1, 0, 0]),
                                ([1, 1, 0])
                             ]
                             )
    def test_spacing_correct(self, spacing):
        mb.Lattice(lattice_spacing=spacing)

    @pytest.mark.parametrize("vectors",
                             [
                                ([[1, 2], [0, 1, 0], [0, 0, 1]]),
                                ([[1, 0, 0], [0, 1, 0], [0, 1, 0]]),
                                (np.identity(4, dtype=np.float64)),
                                ([[1, 2, 3], [3, 2, 1], [2, 1, 3]])
                             ])
    def test_incorrect_lattice_vectors(self, vectors):
        with pytest.raises(ValueError):
            mb.Lattice(lattice_spacing=[1, 1, 1], lattice_vectors=vectors)

    @pytest.mark.parametrize("vectors",
                             [
                                ([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                                ([[1, 0, 0], [-0.5, 0.85, 0], [0, 0, 1]])
                             ])
    def test_correct_lattice_vectors(self, vectors):
        mb.Lattice(lattice_spacing=[1, 1, 1], lattice_vectors=vectors)

    def test_overdefinied_inputs(self):
        space = [1, 1, 1]
        vectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        angles = [90, 90, 90]
        with pytest.raises(ValueError):
            mb.Lattice(lattice_spacing=space, lattice_vectors=vectors,
                       angles=angles)

    @pytest.mark.parametrize("type",
                             [
                                 ([1, 1, 1], 3, list()),
                                 ([1, 1, 1], 3, tuple()),
                                 ([1, 1, 1], 3, str())
                             ]
                             )
    def test_basis_atoms_input_type(self, type):
        with pytest.raises(TypeError):
            mb.Lattice(lattice_spacing=[1, 1, 1], basis_atoms=type)
    # @pytest.mark.parametrize("spacing, dim, x, y, z",
    #                            [
    #                                 ([1, 1, 1], 3, 0, 1, 1),
    #                                 ([1, 1, 1], 3, 1, 0, 1),
    #                                 ([1, 1, 1], 3, 1, 1, 0),
    #                                 ([1, 1, 1], 3, None, 1, 0),
    #                                 ([1, 1, 1], 3, -1, -1, -1)
    #                            ]
    #                          )
    # def test_populate_3d_incorrect_inputs(self, spacing, dim, x, y, z):
    #     test_lattice = mb.Lattice(spacing, dimension=dim)
    #     with pytest.raises(ValueError):
    #         test_lattice.populate(x=x, y=y, z=z)
    #
    # @pytest.mark.parametrize("spacing, dim, x, y",
    #                          [
    #                                 ([1, 1], 2, 0, 1),
    #                                 ([1, 1], 2, 1, 0),
    #                                 ([1, 1], 2, None, 0),
    #                                 ([1, 1], 2, -1, -1)
    #                          ]
    #                          )
    # def test_populate_2d_incorrect_inputs(self, spacing, dim, x, y):
    #     test_lattice = mb.Lattice(spacing, dimension=dim)
    #     with pytest.raises(ValueError):
    #         test_lattice.populate(x=x, y=y)
    #
    # @pytest.mark.parametrize("spacing, dim, x",
    #                          [
    #                                 ([1], 1, 0),
    #                                 ([1], 1, -1)
    #                          ]
    #                          )
    # def test_populate_1d_incorrect_inputs(self, spacing, dim, x):
    #     test_lattice = mb.Lattice(spacing, dimension=dim)
    #     with pytest.raises(ValueError):
    #         test_lattice.populate(x=x)
    #         test_lattice.populate(y=x)
    #         test_lattice.populate(z=x)
    #
    # def test_populate_3d_default(self):
    #     test_lattice = mb.Lattice([1, 1, 1], dimension=3)
    #     a = test_lattice.populate()
    #     np.testing.assert_array_equal(a.xyz[0], [0., 0., 0.])
    #
    # def test_populate_2d_default(self):
    #     test_lattice = mb.Lattice([1, 1], dimension=2)
    #     a = test_lattice.populate()
    #     np.testing.assert_array_equal(a.xyz[0], [0., 0., 0.])
    #
    # def test_populate_1d_default(self):
    #     test_lattice = mb.Lattice([1], dimension=1)
    #     a = test_lattice.populate()
    #     np.testing.assert_array_equal(a.xyz[0], [0., 0., 0.])
    #
    # @pytest.mark.parametrize("my_type",
    #                          [
    #                                 ([]),
    #                                 (()),
    #                                 (np.array),
    #                                 (np.ndarray)
    #                          ]
    #                          )
    # def test_populate_basis_type_incorrect(self, my_type):
    #     test_lattice = mb.Lattice([1, 1, 1], dimension=3)
    #     with pytest.raises(TypeError):
    #         test_lattice.populate(compound_dict=my_type)
    #
    # @pytest.mark.parametrize("not_compound",
    #                          [
    #                                 (1),
    #                                 (mb.Box(lengths=[1, 1, 1])),
    #                                 ("aLattice")
    #                          ]
    #                          )
    # def test_populate_not_compound(self, not_compound):
    #     test_lattice = mb.Lattice([1, 1, 1])
    #     particle_dict = {'default': not_compound}
    #     with pytest.raises(TypeError):
    #         test_lattice.populate(compound_dict=particle_dict)
