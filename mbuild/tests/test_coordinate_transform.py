import numpy as np
import pytest

from mbuild.coordinate_transform import (Translation, CoordinateTransform,
                                         RotationAroundZ, RotationAroundY,
                                         RotationAroundX, Rotation,
                                         ChangeOfBasis, AxisTransform,
                                         RigidTransform, rotate_around_x,
                                         rotate_around_y, rotate_around_z,
                                         force_overlap, translate,
                                         translate_to, x_axis_transform,
                                         y_axis_transform, z_axis_transform,
                                         rotate, spin)
from mbuild.tests.base_test import BaseTest
import mbuild as mb


class TestCoordinateTransform(BaseTest):

    def test_apply_to(self):
        double = CoordinateTransform(T=np.eye(4)*2)
        A = np.array([[1,2,3],[4,5,6],[7,8,9]])
        assert (double.apply_to(A) == np.array([[2,4,6],[8,10,12],[14,16,18]])).all()

    def test_translation(self):
        translation = Translation((10,10,10))
        assert (translation.apply_to(np.array([[1,1,1]])) == np.array([11,11,11])).all()

    def test_rotation_around_z(self):
        Z_rotation = RotationAroundZ(np.pi)
        a = Z_rotation.apply_to(np.array([[2,3,4]]))
        b = np.array([[-2.,-3.,4.]])
        assert np.allclose(a, b, atol=1.e-16)

    def test_rotation_around_y(self):
        Y_rotation = RotationAroundY(np.pi)
        a = Y_rotation.apply_to(np.array([[2,3,4]]))
        b = np.array([-2,3,-4])
        assert np.allclose(a, b, atol=1.e-16)

    def test_rotation_around_x(self):
        X_rotation = RotationAroundX(np.pi)
        a = X_rotation.apply_to(np.array([[2,3,4]]))
        b = np.array([2,-3,-4])
        assert np.allclose(a, b, atol=1.e-16)

    def test_rotation(self):
        rotation = Rotation(np.pi*2/3, np.array([1,1,1]))
        a = rotation.apply_to(np.array([[2,3,4]]))
        b = np.array([4,2,3])
        assert np.allclose(a, b, atol=1.e-16)

    def test_change_of_basis(self):
        change_basis = ChangeOfBasis(np.array([[-2,0,0],[0,-2,0],[0,0,-2]]))
        assert (change_basis.apply_to(np.array([[2,3,4]])) == np.array([[-1.,-1.5,-2.]])).all()

    def test_axis_transform(self):
        origin_transform = AxisTransform(new_origin=np.array([1,1,1]))
        assert (origin_transform.apply_to(np.array([[1,1,1]])) == np.array([[0,0,0]])).all()
        orientation_transform = AxisTransform(point_on_x_axis=np.array([0,0,1]),point_on_xy_plane=np.array([0,1,1]))
        assert (orientation_transform.apply_to(np.array([[2,3,4]])) == np.array([[4,3,-2]])).all()
        axis_transform = AxisTransform(np.array([1,1,1]),np.array([1,1,2]),np.array([1,2,1]))
        assert (axis_transform.apply_to(np.array([[2,3,4]])) == np.array([3,2,-1])).all()

    def test_rigid_transform(self):
        A = np.array([[2,3,4]])
        B = np.array([[3,-4,9]])
        rigid_transform = RigidTransform(A, B)
        assert (rigid_transform.apply_to(np.array([[2,3,4]])) == B).all()

    def test_rotate_0(self, methane):
        before = methane.xyz_with_ports
        rotate(methane, 0.0, np.asarray([1.0, 0.0, 0.0]))
        after = methane.xyz_with_ports
        assert (np.array_equal(before, after))

    def test_rotate_2pi(self, methane):
        before = methane.xyz_with_ports
        rotate(methane, 2*np.pi, np.asarray([1.0, 0.0, 0.0]))
        after = methane.xyz_with_ports
        assert (np.allclose(before, after))

    def test_rotate_zero_vector(self, methane):
        with pytest.raises(ValueError):
            rotate(methane, np.pi/2, np.asarray([0.0, 0.0, 0.0]))
            rotate(methane, 0.0, np.asarray([0.0, 0.0, 0.0]))
            rotate(methane, 2*np.pi, np.asarray([0.0, 0.0, 0.0]))

    def test_spin_zero_vector(self, methane):
        with pytest.raises(ValueError):
            spin(methane, np.pi/2, np.asarray([0.0, 0.0, 0.0]))
            spin(methane, 0.0, np.asarray([0.0, 0.0, 0.0]))
            spin(methane, 2*np.pi, np.asarray([0.0, 0.0, 0.0]))

    def test_spin_360x(self, methane):
        before = methane.xyz_with_ports
        spin(methane, 2*np.pi, np.asarray([1, 0, 0]))
        assert(np.allclose(before, methane.xyz_with_ports, atol=1e-16))

    def test_spin_360y(self, methane):
        before = methane.xyz_with_ports
        spin(methane, 2*np.pi, np.asarray([0, 1, 0]))
        assert(np.allclose(before, methane.xyz_with_ports, atol=1e-16))

    def test_spin_360z(self, methane):
        before = methane.xyz_with_ports
        spin(methane, 2*np.pi, np.asarray([0, 0, 1]))
        assert(np.allclose(before, methane.xyz_with_ports, atol=1e-16))

    def test_spin_0x(self, methane):
        before = methane.xyz_with_ports
        spin(methane, 0, np.asarray([1, 0, 0]))
        assert(np.allclose(before, methane.xyz_with_ports, atol=1e-16))

    def test_spin_0y(self, methane):
        before = methane.xyz_with_ports
        spin(methane, 0, np.asarray([0, 1, 0]))
        assert(np.allclose(before, methane.xyz_with_ports, atol=1e-16))

    def test_spin_0z(self, methane):
        before = methane.xyz_with_ports
        spin(methane, 0, np.asarray([0, 0, 1]))
        assert(np.allclose(before, methane.xyz_with_ports, atol=1e-16))

    def test_spin_x(self, sixpoints):
        before = mb.clone(sixpoints)
        spin(sixpoints, np.pi, np.asarray([1, 0, 0]))
        assert(np.allclose(sixpoints['up'].xyz, before['down'].xyz, atol=1e-16)
                and np.allclose(sixpoints['front'].xyz, before['back'].xyz, atol=1e-16))

    def test_spin_y(self, sixpoints):
        before = mb.clone(sixpoints)
        spin(sixpoints, np.pi, np.asarray([0, 1, 0]))
        assert(np.allclose(sixpoints['left'].xyz, before['right'].xyz, atol=1e-16)
                and np.allclose(sixpoints['front'].xyz, before['back'].xyz, atol=1e-16))

    def test_spin_z(self, sixpoints):
        before = mb.clone(sixpoints)
        spin(sixpoints, np.pi, np.asarray([0, 0, 1]))
        assert(np.allclose(sixpoints['left'].xyz, before['right'].xyz, atol=1e-16)
                and np.allclose(sixpoints['up'].xyz, before['down'].xyz, atol=1e-16))

    def test_spin_arbitraty(self, sixpoints):
        before = mb.clone(sixpoints)
        spin(sixpoints, np.pi, np.asarray([1, 1, 0]))
        assert(np.allclose(sixpoints['up'].xyz, before['right'].xyz, atol=1e-16)
                and np.allclose(sixpoints['down'].xyz, before['left'].xyz, atol=1e-16))

    def test_rotate_around_x(self, methane):
        before = methane.xyz_with_ports
        rotate_around_x(methane, np.pi)
        after = methane.xyz_with_ports
        assert(np.allclose(before[:, 1], -1*after[:, 1], atol=1e-16)
                and np.allclose(before[:, 2], -1*after[:, 2], atol=1e-16))

    def test_rotate_around_y(self, ch2):
        before = ch2.xyz_with_ports
        rotate_around_y(ch2, np.pi)
        after = ch2.xyz_with_ports
        assert(np.allclose(before[:, 0], -1*after[:, 0], atol=1e-16)
                and np.allclose(before[:, 2], -1*after[:, 2], atol=1e-16))

    def test_rotate_around_z(self, ch2):
        before = ch2.xyz_with_ports
        rotate_around_z(ch2, np.pi)
        after = ch2.xyz_with_ports
        assert(np.allclose(before[:, 0], -1*after[:, 0], atol=1e-16)
                and np.allclose(before[:, 1], -1*after[:, 1], atol=1e-16))

    def test_equivalence_transform(self, ch2, ch3, methane):
        ch2_atoms = list(ch2.particles())
        methane_atoms = list(methane.particles())
        force_overlap(ch2, ch2_atoms[0], methane_atoms[0], add_bond=False)
        assert (ch2_atoms[0].pos == methane_atoms[0].pos).all()
        force_overlap(ch2, ch2['up'], ch3['up'])
        assert ch2.n_bonds == 3

        assert ch2.root.bond_graph.number_of_edges() == 3
        assert ch3.root.bond_graph.number_of_edges() == 4

        ethyl = mb.Compound([ch2, ch3])
        assert ethyl.n_bonds == 6

    def test_translate(self, methane):
        methane_atoms = list(methane.particles())
        translate(methane, -methane_atoms[0].pos)
        assert (methane_atoms[0].pos == np.array([0, 0, 0])).all()

    def test_translate_to(self, methane):
        before = methane.xyz_with_ports
        original_center = methane.center
        translate_value = np.array([2, 3, 4])
        translate_to(methane, translate_value)
        assert (methane.xyz_with_ports == before-original_center+translate_value).all()
