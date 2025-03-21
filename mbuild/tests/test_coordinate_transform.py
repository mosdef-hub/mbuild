import itertools as itt

import numpy as np
import pytest

import mbuild as mb
from mbuild.coordinate_transform import (
    AxisTransform,
    ChangeOfBasis,
    CoordinateTransform,
    RigidTransform,
    Rotation,
    RotationAroundX,
    RotationAroundY,
    RotationAroundZ,
    Translation,
    _spin,
    angle,
    force_overlap,
    x_axis_transform,
    y_axis_transform,
    z_axis_transform,
)
from mbuild.tests.base_test import BaseTest


class TestCoordinateTransform(BaseTest):
    def test_apply_to(self):
        double = CoordinateTransform(T=np.eye(4) * 2)
        A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert (
            double.apply_to(A) == np.array([[2, 4, 6], [8, 10, 12], [14, 16, 18]])
        ).all()

    def test_translation(self):
        translation = Translation((10, 10, 10))
        assert (
            translation.apply_to(np.array([[1, 1, 1]])) == np.array([11, 11, 11])
        ).all()

    def test_rotation_around_z(self):
        Z_rotation = RotationAroundZ(np.pi)
        a = Z_rotation.apply_to(np.array([[2, 3, 4]]))
        b = np.array([[-2.0, -3.0, 4.0]])
        assert np.allclose(a, b, atol=1.0e-16)

    def test_rotation_around_y(self):
        Y_rotation = RotationAroundY(np.pi)
        a = Y_rotation.apply_to(np.array([[2, 3, 4]]))
        b = np.array([-2, 3, -4])
        assert np.allclose(a, b, atol=1.0e-16)

    def test_rotation_around_x(self):
        X_rotation = RotationAroundX(np.pi)
        a = X_rotation.apply_to(np.array([[2, 3, 4]]))
        b = np.array([2, -3, -4])
        assert np.allclose(a, b, atol=1.0e-16)

    def test_rotation(self):
        rotation = Rotation(np.pi * 2 / 3, np.array([1, 1, 1]))
        a = rotation.apply_to(np.array([[2, 3, 4]]))
        b = np.array([4, 2, 3])
        assert np.allclose(a, b, atol=1.0e-16)

    def test_change_of_basis(self):
        change_basis = ChangeOfBasis(np.array([[-2, 0, 0], [0, -2, 0], [0, 0, -2]]))
        assert (
            change_basis.apply_to(np.array([[2, 3, 4]]))
            == np.array([[-1.0, -1.5, -2.0]])
        ).all()

    def test_axis_transform(self):
        origin_transform = AxisTransform(new_origin=np.array([1, 1, 1]))
        assert (
            origin_transform.apply_to(np.array([[1, 1, 1]])) == np.array([[0, 0, 0]])
        ).all()
        orientation_transform = AxisTransform(
            point_on_x_axis=np.array([0, 0, 1]),
            point_on_xy_plane=np.array([0, 1, 1]),
        )
        assert (
            orientation_transform.apply_to(np.array([[2, 3, 4]]))
            == np.array([[4, 3, -2]])
        ).all()
        axis_transform = AxisTransform(
            np.array([1, 1, 1]), np.array([1, 1, 2]), np.array([1, 2, 1])
        )
        assert (
            axis_transform.apply_to(np.array([[2, 3, 4]])) == np.array([3, 2, -1])
        ).all()

    def test_rigid_transform(self):
        A = np.array([[2, 3, 4]])
        B = np.array([[3, -4, 9]])
        rigid_transform = RigidTransform(A, B)
        assert (rigid_transform.apply_to(np.array([[2, 3, 4]])) == B).all()

    def test_rotate_0(self, methane):
        before = methane.xyz_with_ports
        methane.rotate(0.0, np.asarray([1.0, 0.0, 0.0]))
        after = methane.xyz_with_ports
        assert np.array_equal(before, after)

    def test_rotate_2pi(self, methane):
        before = methane.xyz_with_ports
        methane.rotate(2 * np.pi, np.asarray([1.0, 0.0, 0.0]))
        after = methane.xyz_with_ports
        assert np.allclose(before, after)

    def test_rotate_zero_vector(self, methane):
        with pytest.raises(ValueError):
            methane.rotate(np.pi / 2, np.asarray([0.0, 0.0, 0.0]))

    def test_spin_zero_vector(self, methane):
        with pytest.raises(ValueError):
            methane.spin(np.pi / 2, np.asarray([0.0, 0.0, 0.0]))

    def test_spin_inputs(self, methane):
        methane.spin(6.9, [1, 0, 0])
        methane.spin(6.9, (1, 0, 0))

    def test_spin_with_anchor(self, methane):
        original_posH = methane[1].pos
        original_posC = methane[0].pos
        methane.spin(6.9, [1, 0, 0], anchor=methane[1])

        assert all(methane[1].pos == original_posH)
        assert any(methane[0].pos != original_posC)

    def test_rotate_inputs(self, methane):
        methane.rotate(6.9, [1, 0, 0])
        methane.rotate(6.9, (1, 0, 0))

    def test_spin_too_many_dimensions_list(self, methane):
        with pytest.raises(ValueError):
            methane.spin(0.1, [1, 0, 0, 0])

    def test_spin_too_many_dimensions_tuple(self, methane):
        with pytest.raises(ValueError):
            methane.spin(0.1, (1, 0, 0, 0))

    def test_rotate_too_many_dimensions_list(self, methane):
        with pytest.raises(ValueError):
            methane.rotate(0.1, [1, 0, 0, 0])

    def test_rotate_too_many_dimensions_tuple(self, methane):
        with pytest.raises(ValueError):
            methane.rotate(0.1, (1, 0, 0, 0))

    def test_spin_too_few_dimensions_list(self, methane):
        with pytest.raises(ValueError):
            methane.spin(0.1, [1, 0])

    def test_spin_too_few_dimensions_tuple(self, methane):
        with pytest.raises(ValueError):
            methane.spin(0.1, (1, 0))

    def test_rotate_too_few_dimensions_list(self, methane):
        with pytest.raises(ValueError):
            methane.rotate(0.1, [1, 0])

    def test_rotate_too_few_dimensions_tuple(self, methane):
        with pytest.raises(ValueError):
            methane.rotate(0.1, (1, 0))

    def test_spin_360x(self, methane):
        before = methane.xyz_with_ports
        methane.spin(2 * np.pi, np.asarray([1, 0, 0]))
        assert np.allclose(before, methane.xyz_with_ports, atol=1e-16)

    def test_spin_360y(self, methane):
        before = methane.xyz_with_ports
        methane.spin(2 * np.pi, np.asarray([0, 1, 0]))
        assert np.allclose(before, methane.xyz_with_ports, atol=1e-16)

    def test_spin_360z(self, methane):
        before = methane.xyz_with_ports
        methane.spin(2 * np.pi, np.asarray([0, 0, 1]))
        assert np.allclose(before, methane.xyz_with_ports, atol=1e-16)

    def test_spin_0x(self, methane):
        before = methane.xyz_with_ports
        methane.spin(0, np.asarray([1, 0, 0]))
        assert np.allclose(before, methane.xyz_with_ports, atol=1e-16)

    def test_spin_0y(self, methane):
        before = methane.xyz_with_ports
        methane.spin(0, np.asarray([0, 1, 0]))
        assert np.allclose(before, methane.xyz_with_ports, atol=1e-16)

    def test_spin_0z(self, methane):
        before = methane.xyz_with_ports
        methane.spin(0, np.asarray([0, 0, 1]))
        assert np.allclose(before, methane.xyz_with_ports, atol=1e-16)

    def test_spin_x(self, sixpoints):
        before = mb.clone(sixpoints)
        sixpoints.spin(np.pi, np.asarray([1, 0, 0]))
        assert np.allclose(
            sixpoints["up"].xyz, before["down"].xyz, atol=1e-16
        ) and np.allclose(sixpoints["front"].xyz, before["back"].xyz, atol=1e-16)

    def test_spin_y(self, sixpoints):
        before = mb.clone(sixpoints)
        sixpoints.spin(np.pi, np.asarray([0, 1, 0]))
        assert np.allclose(
            sixpoints["left"].xyz, before["right"].xyz, atol=1e-16
        ) and np.allclose(sixpoints["front"].xyz, before["back"].xyz, atol=1e-16)

    def test_spin_z(self, sixpoints):
        before = mb.clone(sixpoints)
        sixpoints.spin(np.pi, np.asarray([0, 0, 1]))
        assert np.allclose(
            sixpoints["left"].xyz, before["right"].xyz, atol=1e-16
        ) and np.allclose(sixpoints["up"].xyz, before["down"].xyz, atol=1e-16)

    def test_spin_x_eq(self, sixpoints):
        compound2 = mb.clone(sixpoints)
        sixpoints.spin(np.pi * 1.23456789, np.asarray([1.0, 0.0, 0.0]))

        compound2.spin(np.pi * 1.23456789, around=np.asarray([1, 0, 0]))
        assert np.allclose(compound2.xyz, sixpoints.xyz, atol=1e-16)

    def test_spin_y_eq(self, sixpoints):
        compound2 = mb.clone(sixpoints)
        sixpoints.spin(np.pi * 1.23456789, np.asarray([0.0, 1.0, 0.0]))

        compound2.spin(np.pi * 1.23456789, around=np.asarray([0, 1, 0]))
        assert np.allclose(compound2.xyz, sixpoints.xyz, atol=1e-16)

    def test_spin_z_eq(self, sixpoints):
        compound2 = mb.clone(sixpoints)
        sixpoints.spin(np.pi * 1.23456789, np.asarray([0.0, 0.0, 1.0]))

        compound2.spin(np.pi * 1.23456789, around=np.asarray([0, 0, 1]))
        assert np.allclose(compound2.xyz, sixpoints.xyz, atol=1e-16)

    def test_spin_arbitraty(self, sixpoints):
        before = mb.clone(sixpoints)
        sixpoints.spin(np.pi, np.asarray([1, 1, 0]))
        assert np.allclose(
            sixpoints["up"].xyz, before["right"].xyz, atol=1e-16
        ) and np.allclose(sixpoints["down"].xyz, before["left"].xyz, atol=1e-16)

    def test_spin_relative_compound_coordinates(self, sixpoints):
        """Check compounds's relative coordinates don't change upon spinning"""
        np.random.seed(0)
        angles_before = np.asarray(
            [angle(a, b, c) for (a, b, c) in itt.combinations(sixpoints.xyz, 3)]
        )
        sixpoints.spin(np.pi * 0.1234569789, np.random.rand(3))
        angles_after = np.asarray(
            [angle(a, b, c) for (a, b, c) in itt.combinations(sixpoints.xyz, 3)]
        )
        assert np.allclose(angles_before, angles_after, atol=1e-15)

    def test_equivalence_transform_deprectation_warning(self, ch2):
        ch22 = mb.clone(ch2)
        with pytest.warns(DeprecationWarning):
            mb.equivalence_transform(
                ch22, from_positions=ch22["up"], to_positions=ch2["down"]
            )

    def test_rotate_around_x(self, methane):
        before = methane.xyz_with_ports
        methane.rotate(np.pi, around=np.asarray([1, 0, 0]))
        after = methane.xyz_with_ports
        assert np.allclose(before[:, 1], -1 * after[:, 1], atol=1e-16) and np.allclose(
            before[:, 2], -1 * after[:, 2], atol=1e-16
        )

    def test_rotate_around_y(self, ch2):
        before = ch2.xyz_with_ports
        ch2.rotate(np.pi, around=np.asarray([0, 1, 0]))
        after = ch2.xyz_with_ports
        assert np.allclose(before[:, 0], -1 * after[:, 0], atol=1e-16) and np.allclose(
            before[:, 2], -1 * after[:, 2], atol=1e-16
        )

    def test_rotate_around_z(self, ch2):
        before = ch2.xyz_with_ports
        ch2.rotate(np.pi, around=np.asarray([0, 0, 1]))
        after = ch2.xyz_with_ports
        assert np.allclose(before[:, 0], -1 * after[:, 0], atol=1e-16) and np.allclose(
            before[:, 1], -1 * after[:, 1], atol=1e-16
        )

    def test_rotate_around_x_away_from_origin(self, sixpoints):
        before = sixpoints.xyz_with_ports
        sixpoints.rotate(np.pi, around=np.asarray([1, 0, 0]))
        after = sixpoints.xyz_with_ports
        assert np.allclose(before[:, 1], -1 * after[:, 1], atol=1e-16) and np.allclose(
            before[:, 2], -1 * after[:, 2], atol=1e-16
        )

    def test_rotate_around_y_away_from_origin(self, sixpoints):
        before = sixpoints.xyz_with_ports
        sixpoints.rotate(np.pi, around=np.asarray([0, 1, 0]))
        after = sixpoints.xyz_with_ports
        assert np.allclose(before[:, 0], -1 * after[:, 0], atol=1e-16) and np.allclose(
            before[:, 2], -1 * after[:, 2], atol=1e-16
        )

    def test_rotate_around_z_away_from_origin(self, sixpoints):
        before = sixpoints.xyz_with_ports
        sixpoints.rotate(np.pi, around=np.asarray([0, 0, 1]))
        after = sixpoints.xyz_with_ports
        assert np.allclose(before[:, 1], -1 * after[:, 1], atol=1e-16) and np.allclose(
            before[:, 0], -1 * after[:, 0], atol=1e-16
        )

    def test_rotate_dihedral(self, ethane):
        bond = (ethane[0], ethane[4])
        rotate_angle = np.deg2rad(60)
        ethane.rotate_dihedral(bond, rotate_angle)

        CH_vec1 = ethane[1].pos - ethane[0].pos
        CH_vec2 = ethane[5].pos - ethane[4].pos
        cos_dihedral = np.dot(CH_vec1, CH_vec2) / (
            np.linalg.norm(CH_vec1) * np.linalg.norm(CH_vec2)
        )
        dihedral = np.rad2deg(np.arccos(cos_dihedral))
        assert np.allclose(dihedral, 120, atol=1e-15)

        # Extra test on asymmetric molecule
        compound = mb.load("FCCO", smiles=True)

        # Making sure this is the dihedral we are expecting
        assert compound[0].element.symbol == "F"
        assert compound[1].element.symbol == "C"
        assert compound[2].element.symbol == "C"
        assert compound[3].element.symbol == "O"

        original_dihedral = [
            compound[0].pos,
            compound[1].pos,
            compound[2].pos,
            compound[3].pos,
        ]
        original_angle = dihedral_angle(original_dihedral)

        bond = (compound[1], compound[2])
        rotate_angle = np.deg2rad(68)
        compound.rotate_dihedral(bond, rotate_angle)

        new_dihedral = [
            compound[0].pos,
            compound[1].pos,
            compound[2].pos,
            compound[3].pos,
        ]
        new_angle = dihedral_angle(new_dihedral)

        assert np.allclose(new_angle - original_angle, 68, atol=1e-15)

    def test_equivalence_transform(self, ch2, ch3, methane):
        ch2_atoms = list(ch2.particles())
        methane_atoms = list(methane.particles())
        force_overlap(ch2, ch2_atoms[0], methane_atoms[0], add_bond=False)
        assert (ch2_atoms[0].pos == methane_atoms[0].pos).all()
        force_overlap(ch2, ch2["up"], ch3["up"])
        assert ch2.n_bonds == 3

        assert ch2.root.bond_graph.number_of_edges() == 3
        assert ch3.root.bond_graph.number_of_edges() == 4

        ethyl = mb.Compound([ch2, ch3])
        assert ethyl.n_bonds == 6

    def test_translate(self, methane):
        methane_atoms = list(methane.particles())
        methane.translate(-methane_atoms[0].pos)
        assert (methane_atoms[0].pos == np.array([0, 0, 0])).all()

    def test_translate_to(self, methane):
        before = methane.xyz_with_ports
        original_center = methane.center
        translate_value = np.array([2, 3, 4])
        methane.translate_to(translate_value)
        assert (
            methane.xyz_with_ports == before - original_center + translate_value
        ).all()

    def test_different_translates(self, methane):
        shifted = mb.clone(methane)
        shifted.translate([5, 4, 3])
        shifted_methane_coords = mb.coordinate_transform._translate(
            methane.xyz, [5, 4, 3]
        )
        assert np.array_equal(shifted_methane_coords, shifted.xyz)

    def test_different_translate_tos_origin(self, methane):
        shifted = mb.clone(methane)
        shifted.translate_to([0, 0, 0])
        x = mb.coordinate_transform._translate_to(methane.xyz, [0, 0, 0])
        assert np.array_equal(shifted.xyz, x)

    def test_different_translate_tos_not_origin(self, methane):
        shifted = mb.clone(methane)
        np.random.seed(0)
        point = np.random.rand(3)
        shifted.translate_to(point)
        x = mb.coordinate_transform._translate_to(methane.xyz, point)
        assert np.array_equal(shifted.xyz, x)

    def test_spin(self):
        points = np.asarray(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
            dtype=float,
        )
        new_points_should_be = np.asarray(
            [[0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 0, 0]],
            dtype=float,
        )
        spun_points = _spin(points, np.pi / 2, [0, 0, 1])
        assert np.allclose(spun_points, new_points_should_be, atol=1e-15)

    def test_spin_away_from_origin(self):
        points = np.asarray(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
            dtype=float,
        )
        points += [2, 2, 69]
        new_points_should_be = np.asarray(
            [[2, 2, 69], [2, 3, 69], [1, 2, 69], [2, 1, 69], [3, 2, 69]],
            dtype=float,
        )
        spun_points = _spin(points, np.pi / 2, [0, 0, 1])
        assert np.allclose(spun_points, new_points_should_be, atol=1e-15)

    def test_x_axis_transform(self):
        rot_by_compound = mb.Compound(name="rot_by_compound")
        b = mb.Compound(name="b")
        c = mb.Compound(name="c")
        d = mb.Compound(name="d")
        rot_by_array = mb.Compound(name="rot_by_array")
        b.pos = np.array([0, 0, 0])
        c.pos = np.array([0.5, 0.5, 0.5])
        d.pos = np.array([1, 0, 1])
        array1 = np.array([0, 0, 0])
        array2 = np.array([0.5, 0.5, 0.5])
        array3 = np.array([1, 0, 1])
        x_axis_transform(rot_by_compound, b, c, d)
        x_axis_transform(rot_by_array, array1, array2, array3)
        assert np.array_equal(rot_by_compound.pos, rot_by_array.pos)

    def test_y_axis_transform(self, h2o):
        init_xyz = h2o.xyz
        for i in range(4):
            y_axis_transform(h2o)
        assert np.allclose(h2o.xyz, init_xyz, atol=1e-4)

    def test_z_axis_transform(self, h2o):
        init_xyz = h2o.xyz
        for i in range(4):
            z_axis_transform(h2o)
        assert np.allclose(h2o.xyz, init_xyz, atol=1e-4)

    def test_bondgraph(self, ch3):
        ch3_2 = mb.clone(ch3)
        mb.force_overlap(ch3_2, ch3_2["up"], ch3["up"])
        ch3.add(ch3_2)
        bgraph = ch3.bond_graph
        for edge0, edge1 in bgraph.edges():
            assert bgraph.has_edge(edge0, edge1)
            assert bgraph.has_edge(edge1, edge0)
        neighbors = {"C": 4, "H": 1}
        for node in bgraph.nodes():
            x = map(lambda node: node.name, bgraph._adj[node])
            assert neighbors[node.name] == len(list(x))


# Code to calculate angle of dihedral (stole from https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python#:~:text=The%20angle%20is%20given%20by,Similar%20to%20your%20dihedral2.)
def dihedral_angle(p):
    b1 = p[2] - p[1]
    b0, b1, b2 = -(p[1] - p[0]), b1 / np.sqrt((b1 * b1).sum()), p[3] - p[2]
    v = b0 - (b0[0] * b1[0] + b0[1] * b1[1] + b0[2] * b1[2]) * b1
    w = b2 - (b2[0] * b1[0] + b2[1] * b1[1] + b2[2] * b1[2]) * b1
    x = v[0] * w[0] + v[1] * w[1] + v[2] * w[2]
    y = (
        (b1[1] * v[2] - b1[2] * v[1]) * w[0]
        + (b1[2] * v[0] - b1[0] * v[2]) * w[1]
        + (b1[0] * v[1] - b1[1] * v[0]) * w[2]
    )
    return 180 * np.arctan2(y, x) / np.pi
