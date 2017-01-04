import numpy as np
import pytest
from lattice import Lattice
from mbuild.tests.base_test import BaseTest
import mbuild as mb


class TestLattice(BaseTest):
    """
    Unit Tests for Lattice class functionality.
    """
    def test_dimension_default(self):
        """
        Test the ability of the Lattice class to parse and clean the dimension
        input.
        """
        # default behavior
        a = Lattice(dimension=None)
        assert a.dimension == 3

    def test_dimension_2D(self):
        # 2D system
        a = Lattice(dimension=2)
        assert a.dimension == 2

        # manual setting of 3D
        a = Lattice(dimension=3.0)
        assert a.dimension == 3

        # must be 2D or 3D
        with pytest.raises(ValueError):
            a = Lattice(dimension=1)
            a = Lattice(dimension=4)

        # must be an integer
        with pytest.raises(TypeError):
            a = Lattice(dimension='3')
            a = Lattice(dimension=([1, 2, 3]))

    def test_lattice_vectors(self):
        """
        Unit tests to ensure proper parsing and implementation of the
        lattice_vectors.
        """

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
