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
        spacing = np.reshape(spacing, (3,))
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
                                ([1, 1, 1, 1]),
                                ([1, 'a']),
                                (None),
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

    @pytest.mark.parametrize("the_type",
                             [
                                 (list()),
                                 (tuple()),
                                 (str()),
                                 ([])
                             ]
                             )
    def test_lattice_points_input_type(self, the_type):
        with pytest.raises(TypeError):
            mb.Lattice(lattice_spacing=[1, 1, 1], lattice_points=the_type)

    @pytest.mark.parametrize("incorrect",
                             [
                                 ({'A' : [[.2, .3, .2, .1]]}),
                                 ({'A' : [[None]]}),
                                 ({'A' : [[.2, .3, None]]}),
                                 ({'A' : [[.2, .3, -.5]]}),
                                 ({'A' : [[.2, .3, 1]]}),
                                 ({'A' : [[.2, .3, .1], [.2, .3, .1]]})
                             ]
                             )
    def test_lattice_points_input_type(self, incorrect):
        with pytest.raises(ValueError):
            mb.Lattice(lattice_spacing=[1, 1, 1], lattice_points=incorrect)

    @pytest.mark.parametrize("angles",
                             [
                                ([150, 150, 150]),
                                ([90, 90, -90]),
                                ([90, 90, 180]),
                                ([90, 90, 0]),
                                ([90, 90, 90, 90]),
                                ([97, 3, 120])
                             ]
                             )
    def test_proper_angles(self, angles):
        with pytest.raises(ValueError):
            mb.Lattice(lattice_spacing=[1, 1, 1], angles=angles)

    @pytest.mark.parametrize("vectors, angles",
                             [
                                ([[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                    [90, 90, 90])
                             ]
                             )
    def test_proper_angles(self, vectors, angles):
        testlattice = mb.Lattice(lattice_spacing=[1, 1, 1],
                                  lattice_vectors=vectors)
        np.testing.assert_allclose(testlattice.angles,
                                   np.asarray(angles, dtype=np.float64),
                                   rtol=1e-05, atol=1e-08, equal_nan=False)

    @pytest.mark.parametrize("x, y, z",
                              [
                                (None, 1, 0),
                                (1, None, 1),
                                (1, 1, None),
                                (-1, 1, 1),
                                (1, -1, 1),
                                (1, 1, -1),
                                (1, 1, np.NaN)
                              ])
    def test_incorrect_populate_inputs(self, x, y, z):
        with pytest.raises(ValueError):
            test_lattice = mb.Lattice(lattice_spacing=[1, 1, 1])
            test_lattice.populate(compound_dict={'id':mb.Compound()},
                                  x=x, y=y, z=z)

    @pytest.mark.parametrize("my_type",
                             [
                                ([]),
                                (()),
                                (np.array),
                                (np.ndarray)
                            ]
                            )
    def test_populate_basis_type_incorrect(self, my_type):
        test_lattice = mb.Lattice(lattice_spacing=[1, 1, 1])
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
        test_lattice = mb.Lattice(lattice_spacing=[1, 1, 1])
        particle_dict = {'id': not_compound}
        with pytest.raises(TypeError):
            test_lattice.populate(compound_dict=particle_dict)

    def test_proper_populate(self):
        values_to_check = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                           [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]]
        test_lattice = mb.Lattice(lattice_spacing=[1, 1, 1],
                                  angles=[90, 90, 90])

        new_compound = test_lattice.populate(x=2, y=2, z=2)

        values_to_check = np.asarray(values_to_check, dtype=np.float64)

        is_true = []
        for pos1 in np.split(values_to_check, 8, axis=0):
            for pos2 in np.split(new_compound.xyz, 8, axis=0):
                if np.allclose(pos1, pos2):
                    is_true.append(True)

        assert len(is_true) == len(values_to_check)

    def test_populate_functionalize_no_dict(self):
        # validate port orientation, location, naming
        simple_lat = mb.Lattice([.2, .2, .2], angles=[94.64, 90.67, 96.32])
        expanded_cell = simple_lat.populate(x=3, y=3, z=3, functionalize=True)
        tot, x, y, z, negx, negy, negz = 0, 0, 0, 0, 0, 0, 0
        for ii in expanded_cell.available_ports():
            tot += 1
            n = ii.access_labels
            assert len(n) == 1
            inc1 = n[0].count("negX")
            inc2 = n[0].count("negY")
            inc3 = n[0].count("negZ")
            x += n[0].count("X") - inc1
            y += n[0].count("Y") - inc2
            z += n[0].count("Z") - inc3
            negx += inc1
            negy += inc2
            negz += inc3
            if (inc1 + inc2 + inc3) == 3:
                assert np.allclose(ii.anchor.pos, np.array([0.0, 0.0, 0.0]))
        assert len(expanded_cell.all_ports()) == tot == 26
        assert x == y == z == negx == negy == negz == 9
        assert len(list(expanded_cell.particles())) == 27

    def test_populate_functionalize_complex_unit(self):
        # validate port locations, orientation, naming, and count
        cscl_lattice = mb.Lattice(lattice_spacing=[.4123, .4123, .4123],
                                  lattice_points={'Cl': [[0., 0., 0.]],
                                                  'Cs': [[.5, .5, .5], [.5, .5, 0]]})
        cscl_dict = {'Cl': mb.Compound(name='Cl'),
                     'Cs': mb.Compound(name='Cs')}
        cscl_compound = cscl_lattice.populate(x=4, y=3, z=3,
                                              compound_dict=cscl_dict,
                                              functionalize=True)
        tot, x, y, z, negx, negy, negz = 0, 0, 0, 0, 0, 0, 0
        for ii in cscl_compound.available_ports():
            tot += 1
            n = ii.access_labels
            assert len(n) == 1
            inc1 = n[0].count("negX")
            inc2 = n[0].count("negY")
            inc3 = n[0].count("negZ")
            x += n[0].count("X") - inc1
            y += n[0].count("Y") - inc2
            z += n[0].count("Z") - inc3
            negx += inc1
            negz += inc3
            if (inc1 + inc2 + inc3) == 3 and ii.anchor.name == "CL":
                assert np.allclose(ii.anchor.pos, np.array([0.0, 0.0, 0.0]))
                assert np.allclose(ii.direction/np.linalg.norm(ii.direction),
                                   np.array([-0.57735027, -0.57735027, -0.57735027]))
        assert negz == 30 == y
        assert 24 == z
        assert 23 == x
        assert tot == 72 == len(cscl_compound.all_ports())
        assert len(list(cscl_compound.particles())) == 108

    def test_populate_functionalize_2D_simple_unit(self):
        pass

    def test_populate_functionalize_2D_complex_unit(self):
        pass

    def test_populate_functionalize_triclinic(self):
        pass

    def test_populate_functionalize_validate_basis(self):
        pass

    @pytest.mark.parametrize("no_skin",
                             [
                                 (True)
                             ])
    def test_populate_cannot_skin(self, no_skin):
        pass

    def test_populate_triclinic_skin(self):
        tric_latti = mb.Lattice(lattice_spacing=[1, 1, 1],
                                lattice_points={"Cu": [[0., 0., 0.], [.5, .5, 0.],
                                                       [.5, 0., .5], [0., .5, .5]]},
                                angles=[94.64, 90.67, 96.32])
        comp_dict = {"Cu": mb.Compound(name="Cu")}
        tric_comp = tric_latti.populate(x=3, y=3, z=3,
                                        compound_dict=comp_dict,
                                        functionalize=True)

    def test_set_periodicity(self):
        lattice = mb.Lattice(lattice_spacing=[1, 1, 1], angles=[90, 90, 90],
                             lattice_points={'A': [[0, 0, 0]]})

        compound_test = lattice.populate(compound_dict={'A': mb.Compound()},
                                         x=2, y=5, z=9)

        replication = [2, 5, 9]
        np.testing.assert_allclose(compound_test.periodicity,
                                   np.asarray([x*y for x, y in zip(replication, lattice.lattice_spacing)]))


